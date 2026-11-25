#!/usr/bin/env bash
# Rebuild V1 AUGUSTUS model from the ND90Pr reference.
# Based on the original V1 training procedure supplied by the user.
# Usage: bash train_bipolaris_v1_model.sh [REF_DIR] [TRAIN_DIR]
set -Eeuo pipefail
export LC_ALL=C

REF_DIR="${1:-/mnt/d/Bipolaris_genomes/Bipolaris/GCF_000338995.1}"
TRAIN_DIR="${2:-$HOME/augustus_work/bipolaris_nd90pr}"
GENOME_NAME="GCF_000338995.1_Cocsa1_genomic.fna"
GFF_NAME="genomic.gff"
MODEL="bipolaris_nd90pr"

GENOME="$REF_DIR/$GENOME_NAME"
ANNOTATION="$REF_DIR/$GFF_NAME"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: Activate the Conda environment first." >&2
  exit 1
fi
export AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH:-$CONDA_PREFIX/config}"

for program in gff2gbSmallDNA.pl randomSplit.pl new_species.pl etraining optimize_augustus.pl augustus; do
  command -v "$program" >/dev/null 2>&1 || {
    echo "ERROR: Required program not found: $program" >&2
    exit 1
  }
done

[[ -s "$GENOME" ]] || { echo "ERROR: Genome FASTA missing: $GENOME" >&2; exit 1; }
[[ -s "$ANNOTATION" ]] || { echo "ERROR: Annotation GFF missing: $ANNOTATION" >&2; exit 1; }

SPECIES_DIR="$AUGUSTUS_CONFIG_PATH/species/$MODEL"
if [[ -e "$SPECIES_DIR" ]]; then
  echo "ERROR: V1 model directory already exists: $SPECIES_DIR" >&2
  echo "Nothing was overwritten. Inspect or back it up before retrying." >&2
  exit 1
fi

mkdir -p "$TRAIN_DIR" "$HOME/augustus_work"

GB="$TRAIN_DIR/ND90Pr.genes.gb"
TRAIN_GB="$GB.train"
TEST_GB="$GB.test"

if [[ ! -s "$GB" ]]; then
  echo "Converting GFF and genome FASTA to GenBank training format..."
  gff2gbSmallDNA.pl "$ANNOTATION" "$GENOME" 1000 "$GB" \
    2>&1 | tee "$TRAIN_DIR/gff2gbSmallDNA.log"
fi

if [[ ! -s "$TRAIN_GB" || ! -s "$TEST_GB" ]]; then
  echo "Splitting training examples..."
  randomSplit.pl "$GB" 1000 2>&1 | tee "$TRAIN_DIR/randomSplit.log"
fi

[[ -s "$TRAIN_GB" && -s "$TEST_GB" ]] || {
  echo "ERROR: Training/test split files were not created." >&2
  exit 1
}

echo "Creating AUGUSTUS species model: $MODEL"
new_species.pl --species="$MODEL" 2>&1 | tee "$TRAIN_DIR/new_species.log"
[[ -d "$SPECIES_DIR" ]] || {
  echo "ERROR: new_species.pl did not create $SPECIES_DIR" >&2
  exit 1
}

echo "Running initial etraining. Log: $TRAIN_DIR/etraining.log"
set +e
etraining --species="$MODEL" "$TRAIN_GB" \
  2>&1 | tee "$TRAIN_DIR/etraining.log"
etrain_status=${PIPESTATUS[0]}
set -e
if [[ "$etrain_status" -ne 0 ]]; then
  echo "ERROR: etraining failed with status $etrain_status" >&2
  exit "$etrain_status"
fi

# Work from the Linux filesystem, not the Windows-mounted /mnt/d directory.
WORK="$HOME/augustus_work/${MODEL}_optimization"
mkdir -p "$WORK"
cp "$TRAIN_GB" "$WORK/ND90Pr.genes.gb.train"

echo "Starting optimization. This may take a long time."
echo "Log: $TRAIN_DIR/optimize.log"
cd "$WORK"

set +e
optimize_augustus.pl --species="$MODEL" "$WORK/ND90Pr.genes.gb.train" \
  2>&1 | tee "$TRAIN_DIR/optimize.log"
opt_status=${PIPESTATUS[0]}
set -e
if [[ "$opt_status" -ne 0 ]]; then
  echo "ERROR: optimization failed with status $opt_status" >&2
  exit "$opt_status"
fi

echo "V1 training and optimization completed."
echo "Model: $MODEL"
echo "Model directory: $SPECIES_DIR"
echo "Training directory: $TRAIN_DIR"
echo "Held-out test set: $TEST_GB"
echo "Next: validate this model on $TEST_GB before target-genome comparison."
