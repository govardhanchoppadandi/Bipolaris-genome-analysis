#!/usr/bin/env bash
# Convert reference annotation to AUGUSTUS training data, train, and optimize.
# Usage:
#   bash train_bipolaris_model.sh REF_DIR TRAIN_DIR [GENOME_FASTA] [ANNOTATION_GFF]
set -uo pipefail
export LC_ALL=C

REF_DIR="${1:?Usage: $0 REF_DIR TRAIN_DIR [GENOME_FASTA] [ANNOTATION_GFF]}"
TRAIN_DIR="${2:?Usage: $0 REF_DIR TRAIN_DIR [GENOME_FASTA] [ANNOTATION_GFF]}"
GENOME_NAME="${3:-GCF_000338995.1_Cocsa1_genomic.fna}"
GFF_NAME="${4:-genomic.gff}"
MODEL="${AUGUSTUS_SPECIES:-bipolaris_nd90pr}"

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
  randomSplit.pl "$GB" 1000
fi

if [[ ! -d "$AUGUSTUS_CONFIG_PATH/species/$MODEL" ]]; then
  echo "Creating AUGUSTUS species model: $MODEL"
  new_species.pl --species="$MODEL"
fi

echo "Running etraining. Log: $TRAIN_DIR/etraining.log"
etraining --species="$MODEL" "$TRAIN_GB" \
  2>&1 | tee "$TRAIN_DIR/etraining.log"
etrain_status=${PIPESTATUS[0]}
if [[ "$etrain_status" -ne 0 ]]; then
  echo "ERROR: etraining failed with status $etrain_status" >&2
  exit "$etrain_status"
fi

# Work from the Linux filesystem, not the Windows-mounted /mnt/d directory.
WORK="$HOME/augustus_work/$MODEL"
mkdir -p "$WORK"
cp "$TRAIN_GB" "$WORK/ND90Pr.genes.gb.train"

echo "Starting optimization. This may take a long time."
echo "Log: $TRAIN_DIR/optimize.log"
cd "$WORK" || exit 1
optimize_augustus.pl --species="$MODEL" "$WORK/ND90Pr.genes.gb.train" \
  2>&1 | tee "$TRAIN_DIR/optimize.log"
opt_status=${PIPESTATUS[0]}
if [[ "$opt_status" -ne 0 ]]; then
  echo "ERROR: optimization failed with status $opt_status" >&2
  exit "$opt_status"
fi

echo "Training and optimization command completed for model: $MODEL"
echo "Next: run scripts/validate_bipolaris_model.sh \"$TRAIN_DIR\""
