#!/usr/bin/env bash
# Build a pooled AUGUSTUS V3 model from all 21 user-requested GCA/GCF reference folders.
# GCA/GCF paired accessions may represent duplicate assembly records; all are included
# as requested. The existing V2 and interim V3 models are never changed.
# Selection is based on protein accessions <-> GFF3 CDS protein_id mapping,
# CDS-level completeness, and translation of gffread-extracted CDS sequences.
set -Eeuo pipefail

trap 'rc=$?; echo "ERROR: command failed at line ${LINENO} (exit ${rc}). Check the log files in: ${WORK:-not-created}" >&2; exit "$rc"' ERR

# ---------------------------- USER SETTINGS ----------------------------
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

BASE="/mnt/d/Bipolaris_genomes/Bipolaris"
MODEL="bipolaris_multi_v3_all21"
MAX_GENES_PER_REF=3000
TEST_GENES=1000
FLANK=1000
CPUS=8
OPTIMIZE_ROUNDS=2

# All 21 user-requested reference folders are included.
# NOTE: six GCF accessions are paired records of GCA accessions; including both
# can overweight those assemblies and is NOT equivalent to 21 independent genomes.
# GCA_000527765.1 is included as requested; therefore GCA_000527765.2 is not
# an independent external validation genome for this model.
REFS=(
  GCA_000338975.1
  GCA_000338995.1
  GCA_000354255.1
  GCA_000523435.1
  GCA_000523455.1
  GCA_000527765.1
  GCA_013416765.1
  GCA_019454015.1
  GCA_028564845.1
  GCA_028751825.1
  GCA_028856445.1
  GCA_028858645.1
  GCA_028858655.1
  GCA_028891585.1
  GCA_055773965.2
  GCF_000338995.1
  GCF_000354255.1
  GCF_000523435.1
  GCF_000523455.1
  GCF_000527765.1
  GCF_028858645.1
)
# ----------------------------------------------------------------------

if [[ "${#REFS[@]}" -ne 21 ]]; then
  echo "ERROR: expected exactly 21 reference accessions; found ${#REFS[@]}" >&2
  exit 1
fi

for cmd in augustus etraining optimize_augustus.pl new_species.pl randomSplit.pl gffread python; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: required command not found: $cmd" >&2; exit 1; }
done
GFF2GB="$(command -v gff2gbSmallDNA.pl || true)"
[[ -n "$GFF2GB" ]] || { echo "ERROR: gff2gbSmallDNA.pl is not in PATH" >&2; exit 1; }
[[ -d "$BASE" ]] || { echo "ERROR: base directory does not exist: $BASE" >&2; exit 1; }

SPECIES_DIR="$AUGUSTUS_CONFIG_PATH/species/$MODEL"
if [[ -e "$SPECIES_DIR" ]]; then
  echo "ERROR: model directory already exists: $SPECIES_DIR" >&2
  echo "Nothing was overwritten. Back up/rename that model directory before rebuilding." >&2
  exit 1
fi

# Timestamped work directory prevents accidental mixing with the earlier failed run.
STAMP="$(date +%Y%m%d_%H%M%S)"
WORK="$HOME/augustus_work/${MODEL}_verified_${STAMP}"
mkdir -p "$WORK/per_reference" "$WORK/logs"
echo "Work directory: $WORK"
echo "AUGUSTUS config: $AUGUSTUS_CONFIG_PATH"
echo "Model name: $MODEL (fresh all-21 model; prior models are preserved)"

cat > "$WORK/prepare_one_reference.py" <<'PY'
#!/usr/bin/env python3
import os, re, sys, random
from urllib.parse import unquote
from collections import defaultdict

# args: GFF3 raw.gtf raw.cds.fa genome.fna protein.faa outdir accession max_genes
GFF, RAW_GTF, RAW_CDS, GENOME, PROTEIN, OUTDIR, PREFIX, MAX_GENES = sys.argv[1:]
MAX_GENES = int(MAX_GENES)
os.makedirs(OUTDIR, exist_ok=True)


def parse_gff_attrs(s):
    d = {}
    for item in s.strip().strip(';').split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k.strip()] = unquote(v.strip())
    return d


def safe(s):
    return re.sub(r'[^A-Za-z0-9_.-]', '_', s)


def read_fasta(path):
    """Return ordered list of (full_header, sequence); do not use description as ID."""
    records = []
    header = None
    chunks = []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, ''.join(chunks).upper()))
    return records


# Standard nuclear genetic code. Unknown/ambiguous codons translate to X.
CODON_TEXT = '''
TTT F TTC F TTA L TTG L TCT S TCC S TCA S TCG S TAT Y TAC Y TAA * TAG * TGT C TGC C TGA * TGG W
CTT L CTC L CTA L CTG L CCT P CCC P CCA P CCG P CAT H CAC H CAA Q CAG Q CGT R CGC R CGA R CGG R
ATT I ATC I ATA I ATG M ACT T ACC T ACA T ACG T AAT N AAC N AAA K AAG K AGT S AGC S AGA R AGG R
GTT V GTC V GTA V GTG V GCT A GCC A GCA A GCG A GAT D GAC D GAA E GAG E GGT G GGC G GGA G GGG G
'''.split()
CODON = dict(zip(CODON_TEXT[0::2], CODON_TEXT[1::2]))


def translate_cds(seq):
    if len(seq) % 3:
        return None
    aa = []
    for i in range(0, len(seq), 3):
        aa.append(CODON.get(seq[i:i+3], 'X'))
    return ''.join(aa)


# 1) Read protein FASTA. The first header token is the accession (e.g. EUC26584.1).
protein_seq = {}
protein_header = {}
partial_protein_ids = set()
for header, seq in read_fasta(PROTEIN):
    pid = header.split()[0]
    if pid in protein_seq:
        raise SystemExit(f'ERROR {PREFIX}: duplicate protein accession in protein FASTA: {pid}')
    protein_seq[pid] = seq
    protein_header[pid] = header
    if re.search(r'\bpartial\b', header, flags=re.I):
        partial_protein_ids.add(pid)

if not protein_seq:
    raise SystemExit(f'ERROR {PREFIX}: protein FASTA is empty: {PROTEIN}')

# 2) Parse GFF3. Use CDS protein_id -> CDS Parent transcript -> mRNA Parent gene.
tx_to_gene = {}
protein_to_txs = defaultdict(set)
cds_partial_txs = set()
cds_bad_phase_txs = set()
cds_gff_len = defaultdict(int)
cds_feature_count = defaultdict(int)

with open(GFF, encoding='utf-8', errors='replace') as fh:
    for line in fh:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        typ = f[2].lower()
        a = parse_gff_attrs(f[8])
        if typ in {'mrna', 'transcript'}:
            tid = a.get('ID') or a.get('transcript_id')
            parent = a.get('Parent') or a.get('gene_id') or tid
            if tid and parent:
                tx_to_gene[tid] = parent.split(',')[0]
        elif typ == 'cds':
            parents = a.get('Parent') or a.get('transcript_id')
            pid = a.get('protein_id')
            if not parents or not pid:
                continue
            try:
                seglen = abs(int(f[4]) - int(f[3])) + 1
            except ValueError:
                continue
            phase = f[7]
            for tid in parents.split(','):
                tid = tid.strip()
                if not tid:
                    continue
                protein_to_txs[pid].add(tid)
                cds_gff_len[tid] += seglen
                cds_feature_count[tid] += 1
                # IMPORTANT: only CDS-level partial/range flags are used here.
                # Gene/mRNA partial flags are deliberately NOT used: in
                # GCA_000523435.1 they are present on all genes, while only
                # 1,228 CDS/protein records are actually labelled partial.
                if (a.get('partial', '').lower() in {'true', 'yes', '1'}
                    or 'start_range' in a or 'end_range' in a):
                    cds_partial_txs.add(tid)
                if phase not in {'0', '1', '2'}:
                    cds_bad_phase_txs.add(tid)

# 3) Read gffread-extracted spliced CDS FASTA. gffread -x headers should begin
# with the original transcript ID; this is verified rather than guessed.
cds_by_tx = {}
for header, seq in read_fasta(RAW_CDS):
    tid = header.split()[0]
    if tid in cds_by_tx:
        raise SystemExit(f'ERROR {PREFIX}: duplicate transcript in gffread CDS FASTA: {tid}')
    cds_by_tx[tid] = seq.replace('U', 'T')

if not cds_by_tx:
    raise SystemExit(f'ERROR {PREFIX}: gffread produced no CDS sequences: {RAW_CDS}')

# 4) Select only proteins that are complete-looking in protein FASTA AND whose
# CDS features are not partial, have valid phases, have a complete codon length,
# and translate exactly to the supplied protein sequence.
valid = []
reason = defaultdict(int)
for pid, pseq in protein_seq.items():
    if pid in partial_protein_ids:
        reason['protein FASTA header marked partial'] += 1
        continue
    if not pseq or pseq[0] != 'M':
        reason['protein does not begin with M'] += 1
        continue
    tids = protein_to_txs.get(pid, set())
    if not tids:
        reason['protein accession absent from GFF3 CDS protein_id'] += 1
        continue
    for tid in tids:
        gene = tx_to_gene.get(tid)
        if not gene:
            reason['CDS transcript has no mRNA-to-gene mapping'] += 1
            continue
        if tid in cds_partial_txs:
            reason['CDS feature marked partial/range'] += 1
            continue
        if tid in cds_bad_phase_txs:
            reason['CDS has invalid/missing phase'] += 1
            continue
        cds = cds_by_tx.get(tid)
        if cds is None:
            reason['transcript absent from gffread CDS FASTA'] += 1
            continue
        if len(cds) < 150:
            reason['CDS shorter than 150 bp'] += 1
            continue
        if len(cds) % 3:
            reason['gffread CDS length not divisible by 3'] += 1
            continue
        translated = translate_cds(cds)
        if translated is None:
            reason['CDS translation failed'] += 1
            continue
        if translated.endswith('*'):
            translated = translated[:-1]
        expected = pseq.rstrip('*')
        if translated != expected:
            reason['CDS translation differs from protein FASTA'] += 1
            continue
        valid.append((gene, tid, pid, len(cds), len(pseq)))

if not valid:
    details = '; '.join(f'{k}={v}' for k, v in sorted(reason.items()))
    raise SystemExit(f'ERROR {PREFIX}: zero validated complete protein/CDS models. Filters: {details}')

# 5) Keep one longest validated transcript per gene, avoiding isoform duplication.
best_by_gene = {}
for gene, tid, pid, cds_len, protein_len in valid:
    old = best_by_gene.get(gene)
    if old is None or cds_len > old[3] or (cds_len == old[3] and tid < old[1]):
        best_by_gene[gene] = (gene, tid, pid, cds_len, protein_len)
selected_rows = list(best_by_gene.values())

# Stable, reference-specific random sampling so every genome contributes equally.
seed = 20261102 + sum(ord(c) for c in PREFIX)
rng = random.Random(seed)
if len(selected_rows) > MAX_GENES:
    selected_rows = rng.sample(selected_rows, MAX_GENES)
selected_rows.sort(key=lambda x: x[1])
selected_tids = {row[1] for row in selected_rows}

# 6) Confirm every selected transcript is represented in gffread GTF.
gtf_tids = set()
with open(RAW_GTF, encoding='utf-8', errors='replace') as fh:
    for line in fh:
        m = re.search(r'transcript_id\s+"([^"]+)"', line)
        if m:
            gtf_tids.add(m.group(1))
missing_gtf = selected_tids - gtf_tids
if missing_gtf:
    examples = ', '.join(sorted(missing_gtf)[:8])
    raise SystemExit(f'ERROR {PREFIX}: {len(missing_gtf)} selected transcript IDs absent from GTF. Examples: {examples}')

# 7) Prefix genome sequence IDs and create an accession-safe combined GTF.
seqmap = {}
pref_fa = os.path.join(OUTDIR, PREFIX + '.prefixed.fa')
with open(GENOME, encoding='utf-8', errors='replace') as src, open(pref_fa, 'w') as dst:
    for line in src:
        if line.startswith('>'):
            old = line[1:].strip().split()[0]
            if old in seqmap:
                raise SystemExit(f'ERROR {PREFIX}: duplicate genome FASTA sequence ID: {old}')
            new = safe(PREFIX + '__' + old)
            seqmap[old] = new
            dst.write('>' + new + '\n')
        else:
            dst.write(line.strip() + '\n')

pref_gtf = os.path.join(OUTDIR, PREFIX + '.prefixed.gtf')
good_file = os.path.join(OUTDIR, PREFIX + '.good_transcripts.txt')
missing_seqids = set()
written_lines = 0
with open(RAW_GTF, encoding='utf-8', errors='replace') as src, open(pref_gtf, 'w') as dst:
    for line in src:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        m = re.search(r'transcript_id\s+"([^"]+)"', f[8])
        if not m or m.group(1) not in selected_tids:
            continue
        old_seqid = f[0]
        if old_seqid not in seqmap:
            missing_seqids.add(old_seqid)
            continue
        f[0] = seqmap[old_seqid]
        attrs = f[8]
        attrs = re.sub(r'(transcript_id\s+")([^"]+)(")',
                       lambda m: m.group(1) + safe(PREFIX + '__' + m.group(2)) + m.group(3), attrs)
        attrs = re.sub(r'(gene_id\s+")([^"]+)(")',
                       lambda m: m.group(1) + safe(PREFIX + '__' + m.group(2)) + m.group(3), attrs)
        f[8] = attrs
        dst.write('\t'.join(f) + '\n')
        written_lines += 1

if missing_seqids:
    raise SystemExit(f'ERROR {PREFIX}: selected GTF records use sequence IDs absent from genome FASTA: {sorted(missing_seqids)[:8]}')
if written_lines == 0:
    raise SystemExit(f'ERROR {PREFIX}: no selected GTF rows were written.')

with open(good_file, 'w') as fh:
    for tid in sorted(selected_tids):
        fh.write(safe(PREFIX + '__' + tid) + '\n')

mapping_file = os.path.join(OUTDIR, PREFIX + '.selected_gene_models.tsv')
with open(mapping_file, 'w') as fh:
    fh.write('gene_id\ttranscript_id\tprotein_accession\tcds_bp\tprotein_aa\n')
    for gene, tid, pid, cds_len, protein_len in selected_rows:
        fh.write(f'{gene}\t{tid}\t{pid}\t{cds_len}\t{protein_len}\n')

print(f'{PREFIX}: proteins={len(protein_seq)}; protein-partial={len(partial_protein_ids)}')
print(f'{PREFIX}: validated complete protein/CDS transcript models={len(valid)}')
print(f'{PREFIX}: unique genes after one-isoform-per-gene={len(best_by_gene)}')
print(f'{PREFIX}: selected for training={len(selected_rows)}')
print(f'{PREFIX}: GTF rows retained={written_lines}')
print(f'{PREFIX}: filtered counts: ' + '; '.join(f'{k}={v}' for k, v in sorted(reason.items())))
print(f'{PREFIX}: selected mapping={mapping_file}')
print(f'{PREFIX}: prefixed genome={pref_fa}')
print(f'{PREFIX}: prefixed GTF={pref_gtf}')
print(f'{PREFIX}: good transcript list={good_file}')
PY

chmod +x "$WORK/prepare_one_reference.py"

ALL_GB="$WORK/${MODEL}.all_references.gb"
: > "$ALL_GB"
TOTAL_RECORDS=0
USED_REFS=0
: > "$WORK/reference_summary.tsv"
printf 'accession\tvalidated_selected\tGenBank_records\tstatus\n' >> "$WORK/reference_summary.tsv"

for ID in "${REFS[@]}"; do
  DIR="$BASE/$ID"
  OUT="$WORK/per_reference/$ID"
  mkdir -p "$OUT"

  # The accession-prefixed pattern prevents accidentally selecting cds_from_genomic.fna.
  GENOME="$(find "$DIR" -maxdepth 1 -type f -name "${ID}*_genomic.fna" -print -quit 2>/dev/null || true)"
  GFF="$(find "$DIR" -maxdepth 1 -type f \( -name 'genomic.gff' -o -name 'genomic.gff3' \) -print -quit 2>/dev/null || true)"
  if [[ -z "$GFF" ]]; then
    GFF="$(find "$DIR" -maxdepth 1 -type f \( -name '*.gff' -o -name '*.gff3' \) -print -quit 2>/dev/null || true)"
  fi
  PROTEIN="$(find "$DIR" -maxdepth 1 -type f \( -name 'protein.faa' -o -name "${ID}*protein*.faa" \) -print -quit 2>/dev/null || true)"

  if [[ -z "$GENOME" || -z "$GFF" || -z "$PROTEIN" || ! -s "$GENOME" || ! -s "$GFF" || ! -s "$PROTEIN" ]]; then
    echo "ERROR: required genome/GFF/protein file missing for $ID" | tee -a "$WORK/logs/missing_inputs.log"
    printf '%s\t0\t0\tMISSING_INPUTS\n' "$ID" >> "$WORK/reference_summary.tsv"
    echo "STOP: all 21 references are required; correct the missing files and rerun." >&2
    exit 1
  fi

  echo
  echo "================ Preparing $ID ================"
  echo "Genome:  $GENOME"
  echo "GFF3:    $GFF"
  echo "Protein: $PROTEIN"

  RAW_GTF="$OUT/$ID.raw.gtf"
  RAW_CDS="$OUT/$ID.raw.cds.fa"
  if ! gffread -E -g "$GENOME" -T -o "$RAW_GTF" -x "$RAW_CDS" "$GFF" > "$OUT/gffread.stdout.log" 2> "$OUT/gffread.stderr.log"; then
    echo "ERROR: gffread failed for $ID; see $OUT/gffread.stderr.log" | tee -a "$WORK/logs/reference_failures.log"
    printf '%s\t0\t0\tGFFREAD_FAILED\n' "$ID" >> "$WORK/reference_summary.tsv"
    exit 1
  fi
  if [[ ! -s "$RAW_GTF" || ! -s "$RAW_CDS" ]]; then
    echo "ERROR: gffread produced empty GTF/CDS for $ID" | tee -a "$WORK/logs/reference_failures.log"
    printf '%s\t0\t0\tEMPTY_GFFREAD_OUTPUT\n' "$ID" >> "$WORK/reference_summary.tsv"
    exit 1
  fi

  if ! python "$WORK/prepare_one_reference.py" "$GFF" "$RAW_GTF" "$RAW_CDS" "$GENOME" "$PROTEIN" "$OUT" "$ID" "$MAX_GENES_PER_REF" 2>&1 | tee "$OUT/preparation.log"; then
    echo "ERROR: preparation failed for $ID; see $OUT/preparation.log" | tee -a "$WORK/logs/reference_failures.log"
    printf '%s\t0\t0\tPREPARATION_FAILED\n' "$ID" >> "$WORK/reference_summary.tsv"
    exit 1
  fi

  PREF_GTF="$OUT/$ID.prefixed.gtf"
  PREF_FA="$OUT/$ID.prefixed.fa"
  GOOD="$OUT/$ID.good_transcripts.txt"
  GB="$OUT/$ID.training.gb"

  if ! "$GFF2GB" --good="$GOOD" "$PREF_GTF" "$PREF_FA" "$FLANK" "$GB" > "$OUT/gff2gb.stdout.log" 2> "$OUT/gff2gb.stderr.log"; then
    echo "ERROR: gff2gbSmallDNA failed for $ID; see $OUT/gff2gb.stderr.log" | tee -a "$WORK/logs/reference_failures.log"
    printf '%s\t0\t0\tGFF2GB_FAILED\n' "$ID" >> "$WORK/reference_summary.tsv"
    exit 1
  fi

  NREC="$(grep -c '^LOCUS' "$GB" || true)"
  NSEL="$(wc -l < "$GOOD" | tr -d ' ')"
  echo "$ID: selected transcript IDs=$NSEL; GenBank records=$NREC"
  if [[ "$NREC" -lt 100 ]]; then
    echo "ERROR: fewer than 100 GenBank loci for $ID; see $OUT/gff2gb.stderr.log" | tee -a "$WORK/logs/reference_failures.log"
    printf '%s\t%s\t%s\tTOO_FEW_GENBANK_RECORDS\n' "$ID" "$NSEL" "$NREC" >> "$WORK/reference_summary.tsv"
    exit 1
  fi

  cat "$GB" >> "$ALL_GB"
  TOTAL_RECORDS=$((TOTAL_RECORDS + NREC))
  USED_REFS=$((USED_REFS + 1))
  printf '%s\t%s\t%s\tUSED\n' "$ID" "$NSEL" "$NREC" >> "$WORK/reference_summary.tsv"
done

printf '\n================ Preparation summary ================\n'
cat "$WORK/reference_summary.tsv"
printf '\nReferences used: %s\nTotal GenBank loci: %s\n' "$USED_REFS" "$TOTAL_RECORDS"

# Require enough independent annotated references and loci to avoid building
# a misleadingly small pooled model if several references failed.
if [[ "$USED_REFS" -ne 21 || "$TOTAL_RECORDS" -lt 5000 ]]; then
  echo "ERROR: expected all 21 references; only $USED_REFS were used ($TOTAL_RECORDS loci). Model was NOT created." >&2
  echo "Review $WORK/reference_summary.tsv and the per-reference logs." >&2
  exit 1
fi

printf 'References used: %s\nTotal training loci: %s\n' "$USED_REFS" "$TOTAL_RECORDS" > "$WORK/training_set_summary.txt"

cd "$WORK"
randomSplit.pl "$ALL_GB" "$TEST_GENES" > "$WORK/randomSplit.log" 2>&1
TRAIN="$ALL_GB.train"
TEST="$ALL_GB.test"
[[ -s "$TRAIN" && -s "$TEST" ]] || { echo "ERROR: randomSplit.pl did not create non-empty train/test files; see $WORK/randomSplit.log" >&2; exit 1; }
NTRAIN="$(grep -c '^LOCUS' "$TRAIN" || true)"
NTEST="$(grep -c '^LOCUS' "$TEST" || true)"
if [[ "$NTRAIN" -lt 1000 || "$NTEST" -lt 100 ]]; then
  echo "ERROR: unexpected train/test split: train=$NTRAIN test=$NTEST" >&2
  exit 1
fi
printf 'Training loci: %s\nTest loci: %s\n' "$NTRAIN" "$NTEST" | tee "$WORK/split_summary.txt"

# Create a new species model. V2 remains untouched.
new_species.pl --species="$MODEL" > "$WORK/new_species.log" 2>&1
[[ -d "$SPECIES_DIR" ]] || { echo "ERROR: new_species.pl did not create $SPECIES_DIR" >&2; exit 1; }

# Initial training.
etraining --species="$MODEL" "$TRAIN" > "$WORK/etraining_initial.log" 2>&1

# Optimize the AUGUSTUS meta-parameters using the training partition only.
optimize_augustus.pl \
  --species="$MODEL" \
  --metapars="$SPECIES_DIR/${MODEL}_metapars.cfg" \
  --rounds="$OPTIMIZE_ROUNDS" \
  --cpus="$CPUS" \
  --cleanup=1 \
  "$TRAIN" > "$WORK/optimize.log" 2>&1

# Final parameter estimation on training partition.
etraining --species="$MODEL" "$TRAIN" > "$WORK/etraining_final.log" 2>&1

# Save a separate backup copy of the finished species directory.
cp -a "$SPECIES_DIR" "$WORK/model_copy"

cat > "$WORK/README_model.txt" <<EOF
AUGUSTUS pooled model: $MODEL
Model directory: $SPECIES_DIR
Training work directory: $WORK
References used: $USED_REFS
Total GenBank loci before split: $TOTAL_RECORDS
Training loci: $NTRAIN
Held-out random loci: $NTEST
Maximum selected genes per reference: $MAX_GENES_PER_REF

All 21 requested accession folders were included. Some GCA/GCF pairs may
represent duplicate assembly records and therefore can overweight those
assemblies; these are not 21 independent biological genomes.
The held-out partition is a random gene-level split pooled across all records.
It is not a leave-one-genome-out or leave-one-species-out validation.
GCA_000527765.1 is included, so GCA_000527765.2 is related to a training
assembly and must not be described as an independent external validation.
EOF

echo
echo "================ V3 MODEL TRAINING FINISHED ================"
echo "Model name:       $MODEL"
echo "Model directory:  $SPECIES_DIR"
echo "Work directory:   $WORK"
echo "Reference folders used: $USED_REFS (all 21 required)"
echo "Training loci:    $NTRAIN"
echo "Held-out loci:    $NTEST"
echo "Reference report: $WORK/reference_summary.tsv"
echo "Optimization log: $WORK/optimize.log"
echo "Final training:   $WORK/etraining_final.log"
echo "No predictions were run on GCA_000527765.2."
