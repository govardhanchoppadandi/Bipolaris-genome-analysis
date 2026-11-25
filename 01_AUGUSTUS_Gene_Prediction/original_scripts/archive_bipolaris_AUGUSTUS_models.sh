
#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# BIPOLARIS AUGUSTUS MODEL PROVENANCE ARCHIVE
# Archives V1-V5 models, available scripts, logs, configurations,
# training-data inventory, checksums, and model-selection report.
# ============================================================

ROOT="/mnt/d/Bipolaris_genomes/Bipolaris"
ARI_ROOT="/mnt/d/Bipolaris_genomes/ARI genomes"

ARCHIVE="$ROOT/AUGUSTUS_V1_V5_COMPLETE_ARCHIVE"
CONDA_CONFIG="${CONDA_PREFIX:-$HOME/miniconda3/envs/busco_latest}/config"

mkdir -p "$ARCHIVE"/{models,code_and_logs,comparison_results,training_data_inventory,metadata}

echo "Starting AUGUSTUS model archive..."
echo "Archive location: $ARCHIVE"

# ------------------------------------------------------------
# 1. Record system and software information
# ------------------------------------------------------------

{
    echo "Bipolaris AUGUSTUS model archive"
    echo "Date: $(date -Is)"
    echo "Computer: $(hostname)"
    echo "Operating system:"
    uname -a
    echo
    echo "Conda environment:"
    conda info --envs 2>/dev/null || true
    echo
    echo "Software versions:"
    augustus --version 2>&1 || true
    busco --version 2>&1 || true
    gffread --version 2>&1 || true
    echo
    echo "AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH:-not_set}"
    echo "Conda config path: $CONDA_CONFIG"
} > "$ARCHIVE/metadata/software_environment.txt"

# ------------------------------------------------------------
# 2. Copy all five trained AUGUSTUS species models
# ------------------------------------------------------------

MODELS=(
    "bipolaris_nd90pr"
    "bipolaris_nd90pr_v2"
    "bipolaris_multi_v3_all21"
    "bipolaris_multi_v4_advanced"
    "bipolaris_multi_v5_clean"
)

for MODEL in "${MODELS[@]}"; do
    SOURCE="$CONDA_CONFIG/species/$MODEL"

    if [[ -d "$SOURCE" ]]; then
        echo "Copying model: $MODEL"
        mkdir -p "$ARCHIVE/models/$MODEL"
        cp -a "$SOURCE/." "$ARCHIVE/models/$MODEL/"
    else
        echo "WARNING: Model directory not found: $SOURCE"
        echo "$MODEL NOT FOUND: $SOURCE" \
            >> "$ARCHIVE/metadata/missing_items.txt"
    fi
done

# Also record the V5 deduplicated alternative if present.
ALT_MODEL="$CONDA_CONFIG/species/bipolaris_multi_v5_deduplicated"
if [[ -d "$ALT_MODEL" ]]; then
    cp -a "$ALT_MODEL" "$ARCHIVE/models/"
fi

# ------------------------------------------------------------
# 3. Find and copy available training scripts, logs and configs
# ------------------------------------------------------------

SEARCH_ROOTS=(
    "$ROOT"
    "$ARI_ROOT"
    "$HOME"
)

for SEARCH in "${SEARCH_ROOTS[@]}"; do
    [[ -d "$SEARCH" ]] || continue

    find "$SEARCH" \
        -maxdepth 6 \
        -type f \
        \( -iname "*.sh" \
        -o -iname "*.log" \
        -o -iname "*.cfg" \
        -o -iname "*.params" \
        -o -iname "*.out" \
        -o -iname "*.err" \
        -o -iname "*.txt" \
        -o -iname "*.yaml" \
        -o -iname "*.yml" \) \
        ! -path "$ARCHIVE/*" \
        ! -path "*/.cache/*" \
        ! -path "*/site-packages/*" \
        ! -path "*/conda/pkgs/*" \
        -print0 |
    while IFS= read -r -d '' FILE; do
        LOWER="$(basename "$FILE" | tr '[:upper:]' '[:lower:]')"

        # Save likely model-development files, not every unrelated file.
        if [[ "$LOWER" =~ august|bipolaris|training|train|etraining|optimi|busco|v[1-5] ]]; then
            cp --parents -n "$FILE" "$ARCHIVE/code_and_logs/" 2>/dev/null || true
        fi
    done
done

# Save shell history if available.
for HIST in "$HOME/.bash_history" "$HOME/.zsh_history"; do
    if [[ -f "$HIST" ]]; then
        cp -n "$HIST" "$ARCHIVE/code_and_logs/" || true
    fi
done

# ------------------------------------------------------------
# 4. Preserve the final V1-V5 comparison results
# ------------------------------------------------------------

COMPARISON="$ROOT/FINAL_V1_V2_V3_COMPARISON"

if [[ -d "$COMPARISON" ]]; then
    cp -a "$COMPARISON" "$ARCHIVE/comparison_results/"
fi

# Preserve V2 batch prediction logs and report if they exist.
for ITEM in \
    "$ROOT/V2_AUGUSTUS_LOGS" \
    "$ROOT/V2_gene_prediction_report.tsv"; do
    if [[ -e "$ITEM" ]]; then
        cp -a "$ITEM" "$ARCHIVE/comparison_results/"
    fi
done

# ------------------------------------------------------------
# 5. Inventory all possible training input data
# ------------------------------------------------------------

INVENTORY="$ARCHIVE/training_data_inventory/training_data_file_inventory.tsv"

printf "file_path\tfile_size_bytes\tsha256\n" > "$INVENTORY"

for SEARCH in "$ROOT" "$ARI_ROOT"; do
    [[ -d "$SEARCH" ]] || continue

    find "$SEARCH" \
        -type f \
        \( -iname "*.fna" \
        -o -iname "*.fasta" \
        -o -iname "*.fa" \
        -o -iname "*.gff" \
        -o -iname "*.gff3" \
        -o -iname "*.gbk" \
        -o -iname "*.gbff" \
        -o -iname "*.genbank" \
        -o -iname "*.faa" \
        -o -iname "*.pep" \
        -o -iname "*.cds" \
        -o -iname "*.gff.gz" \
        -o -iname "*.fna.gz" \
        -o -iname "*.fasta.gz" \) \
        -print0 |
    while IFS= read -r -d '' FILE; do
        SIZE="$(stat -c %s "$FILE")"
        HASH="$(sha256sum "$FILE" | awk '{print $1}')"
        printf '%s\t%s\t%s\n' "$FILE" "$SIZE" "$HASH" \
            >> "$INVENTORY"
    done
done

# ------------------------------------------------------------
# 6. Save the known five-model BUSCO comparison
# ------------------------------------------------------------

cat > "$ARCHIVE/metadata/V1_V5_BUSCO_comparison.tsv" <<'EOF'
Model	Genes	BUSCO_Complete	Single_Copy	Duplicated	Fragmented	Missing
V1_bipolaris_nd90pr	9994	93.7	93.6	0.1	1.1	5.2
V2_bipolaris_nd90pr_v2	10196	94.9	94.7	0.1	1.3	3.8
V3_bipolaris_multi_v3_all21	9931	93.8	93.7	0.1	1.1	5.1
V4_bipolaris_multi_v4_advanced	9684	92.5	92.4	0.1	1.2	6.3
V5_bipolaris_multi_v5_clean	9911	93.8	93.7	0.1	1.1	5.1
EOF

# ------------------------------------------------------------
# 7. Document why V2 was selected
# ------------------------------------------------------------

cat > "$ARCHIVE/metadata/WHY_V2_WAS_SELECTED.txt" <<'EOF'
MODEL SELECTION REPORT: BIPOLARIS AUGUSTUS MODELS V1-V5

Evaluation:
All five AUGUSTUS models were evaluated on the same test genome:
GCA_000527765.2_ANU_Cvic_canu1.6_pilon5_genomic.fna

Protein predictions were evaluated using:
BUSCO 6.1.0
Lineage: pleosporales_odb12
BUSCO groups: 6066

RESULTS:

V1:
Genes: 9994
Complete BUSCO: 93.7%
Single-copy: 93.6%
Duplicated: 0.1%
Fragmented: 1.1%
Missing: 5.2%

V2:
Genes: 10196
Complete BUSCO: 94.9%
Single-copy: 94.7%
Duplicated: 0.1%
Fragmented: 1.3%
Missing: 3.8%

V3:
Genes: 9931
Complete BUSCO: 93.8%
Single-copy: 93.7%
Duplicated: 0.1%
Fragmented: 1.1%
Missing: 5.1%

V4:
Genes: 9684
Complete BUSCO: 92.5%
Single-copy: 92.4%
Duplicated: 0.1%
Fragmented: 1.2%
Missing: 6.3%

V5:
Genes: 9911
Complete BUSCO: 93.8%
Single-copy: 93.7%
Duplicated: 0.1%
Fragmented: 1.1%
Missing: 5.1%

REASON FOR SELECTING V2:

1. V2 had the highest complete BUSCO score: 94.9%.
2. V2 had the highest single-copy BUSCO score: 94.7%.
3. V2 had the lowest missing BUSCO percentage: 3.8%.
4. V2 predicted 10,196 genes.
5. BUSCO duplicated percentage was low (0.1%).

CONCLUSION:

V2 (bipolaris_nd90pr_v2) was selected as the current
best-performing AUGUSTUS model based on protein BUSCO
completeness on the common test genome.

IMPORTANT LIMITATIONS:

- V2 is the best model by the reported BUSCO comparison,
  not necessarily the most accurate model for every gene.
- V2 has a slightly higher fragmented BUSCO percentage
  (1.3%) than V1, V3 and V5 (1.1%).
- Gene structure accuracy, exon boundaries, intron boundaries,
  UTRs and false-positive predictions were not fully assessed
  by this BUSCO comparison.
- Independent validation with trusted annotations or RNA-seq
  evidence is recommended.
- The test genome GCA_000527765.2 is not completely independent
  of the V3 training set because GCA_000527765.1 was included
  in the V3 reference set.
- Whole-genome BUSCO results are assembly-level results and
  must not be directly compared with protein-prediction BUSCO.

Therefore, V2 is the current selected model based on the
available evaluation, subject to further structural validation.
EOF

# ------------------------------------------------------------
# 8. Save model directory contents
# ------------------------------------------------------------

find "$ARCHIVE/models" -type f -printf '%P\t%s bytes\n' \
    > "$ARCHIVE/metadata/model_files_inventory.txt"

# ------------------------------------------------------------
# 9. Generate checksums for archived files
# ------------------------------------------------------------

(
    cd "$ARCHIVE"
    find . -type f ! -name "ARCHIVE_SHA256SUMS.txt" -print0 |
    sort -z |
    xargs -0 -r sha256sum
) > "$ARCHIVE/metadata/ARCHIVE_SHA256SUMS.txt"

# ------------------------------------------------------------
# 10. Create archive README
# ------------------------------------------------------------

cat > "$ARCHIVE/README.txt" <<'EOF'
BIPOLARIS AUGUSTUS V1-V5 MODEL ARCHIVE

Contents:

models/
    Copies of the trained AUGUSTUS species model directories.

code_and_logs/
    Training scripts, optimization scripts, logs, configuration
    files and shell history found during archive creation.

comparison_results/
    Available V1-V5 comparison results and V2 prediction reports.

training_data_inventory/
    Inventory of potential training data files with file sizes
    and SHA-256 checksums.

metadata/
    Software environment information.
    Model BUSCO comparison.
    V2 model selection justification.
    Model file inventory.
    Archive checksums.

IMPORTANT:
The training data inventory identifies files available in the
project folders. It does not prove that every listed file was
used to train each model. Confirm actual training inputs using
the original training scripts and logs.
EOF

# ------------------------------------------------------------
# 11. Finish
# ------------------------------------------------------------

echo
echo "=================================================="
echo "ARCHIVE COMPLETED"
echo "=================================================="
echo "Location: $ARCHIVE"
echo
du -sh "$ARCHIVE"
echo
echo "Files saved:"
find "$ARCHIVE" -maxdepth 2 -type f | sort
echo
echo "To compress the archive, run:"
echo "tar -czf '${ARCHIVE}.tar.gz' -C '$ROOT' 'AUGUSTUS_V1_V5_COMPLETE_ARCHIVE'"
