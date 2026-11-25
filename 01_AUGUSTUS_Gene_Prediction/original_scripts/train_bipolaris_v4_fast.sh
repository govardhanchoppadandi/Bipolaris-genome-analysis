#!/usr/bin/env bash
set -Eeuo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest

export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

SPECIES="bipolaris_multi_v4_advanced"

OLD_WORK="$HOME/augustus_work/bipolaris_multi_v3_all21_verified_20261102_100949"

WORK="$HOME/augustus_work/${SPECIES}_fast"

SPEC_DIR="$AUGUSTUS_CONFIG_PATH/species/$SPECIES"

FULL_TRAIN="$OLD_WORK/bipolaris_multi_v3_all21.all_references.gb.train"

TEST="$OLD_WORK/bipolaris_multi_v3_all21.all_references.gb.test"

mkdir -p "$WORK"

LOG="$WORK/V4_training.log"

exec > >(tee -a "$LOG") 2>&1

echo "============================================"
echo "BIPOLARIS V4 ADVANCED FAST TRAINING"
echo "============================================"

# Safety checks
if [[ ! -s "$FULL_TRAIN" ]]; then
    echo "ERROR: V3 training file not found."
    exit 1
fi

if [[ ! -s "$TEST" ]]; then
    echo "ERROR: V3 test file not found."
    exit 1
fi

if [[ -d "$SPEC_DIR" ]]; then
    echo "ERROR: V4 species directory already exists."
    echo "No existing model will be overwritten."
    exit 1
fi

echo
echo "========== SYSTEM =========="
echo "CPU threads: $(nproc)"
free -h

CPUS=$(nproc)
if (( CPUS > 8 )); then CPUS=8; fi

echo "Optimization CPUs: $CPUS"

echo
echo "========== INPUT TRAINING DATA =========="
ls -lh "$FULL_TRAIN"

# Create a stratified optimization sample:
# 250 randomly selected records from every accession.

echo
echo "========== CREATE OPTIMIZATION SAMPLE =========="

export FULL_TRAIN WORK

python - <<'PY'
from Bio import SeqIO
from collections import defaultdict
import random
import os

source = os.environ["FULL_TRAIN"]
work = os.environ["WORK"]

random.seed(20261110)

counts = defaultdict(list)

print("Scanning accession records...")

for i, rec in enumerate(SeqIO.parse(source, "genbank")):
    accession = rec.id.split("__")[0]
    counts[accession].append(i)

print("Accessions found:", len(counts))

selected = set()

for accession, indices in sorted(counts.items()):
    chosen = random.sample(indices, min(250, len(indices)))
    selected.update(chosen)
    print(accession, "selected:", len(chosen))

sample = os.path.join(work, "V4_optimization_sample.gb")

print("\nWriting optimization sample...")

written = 0

with open(sample, "w") as out:
    for i, rec in enumerate(SeqIO.parse(source, "genbank")):
        if i in selected:
            SeqIO.write(rec, out, "genbank")
            written += 1

print("\nTotal selected:", written)
print("Sample:", sample)

if len(counts) != 21:
    raise SystemExit(
        f"STOP: Expected 21 accession groups, found {len(counts)}"
    )
PY

SAMPLE="$WORK/V4_optimization_sample.gb"

if [[ ! -s "$SAMPLE" ]]; then
    echo "ERROR: Optimization sample missing."
    exit 1
fi

echo
echo "========== SAMPLE SIZE =========="
grep -c '^LOCUS' "$SAMPLE"

# Create a completely separate AUGUSTUS species.

echo
echo "========== CREATE V4 SPECIES =========="

new_species.pl --species="$SPECIES"

# Initial training on optimization subset.

echo
echo "========== INITIAL V4 TRAINING =========="

etraining \
    --species="$SPECIES" \
    "$SAMPLE" \
    > "$WORK/etraining_initial.log" 2>&1

# Optimize using the representative subset.

echo
echo "========== V4 OPTIMIZATION: 2 ROUNDS =========="

optimize_augustus.pl \
    --species="$SPECIES" \
    --rounds=2 \
    --cpus="$CPUS" \
    "$SAMPLE" \
    > "$WORK/optimize.log" 2>&1

# Final training using every original training record.

echo
echo "========== FINAL V4 TRAINING: ALL RECORDS =========="

etraining \
    --species="$SPECIES" \
    "$FULL_TRAIN" \
    > "$WORK/etraining_final.log" 2>&1

# Evaluate on the same pooled test set used for V3.

echo
echo "========== V4 TEST EVALUATION =========="

augustus \
    --species="$SPECIES" \
    --UTR=off \
    "$TEST" \
    > "$WORK/V4_test_predictions.out" \
    2> "$WORK/V4_test_evaluation.log"

echo
echo "========== V4 ACCURACY =========="

grep -A35 "Evaluation of gene prediction" \
    "$WORK/V4_test_predictions.out" | head -35

echo
echo "========== V4 COMPLETED =========="

echo "Species: $SPECIES"
echo "Model: $SPEC_DIR"
echo "Training workspace: $WORK"

echo "============================================"
echo "V4 TRAINING AND INITIAL EVALUATION FINISHED"
echo "============================================"
