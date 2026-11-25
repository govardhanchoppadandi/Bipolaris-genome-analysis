#!/usr/bin/env bash
# Validate a trained AUGUSTUS model against its held-out test set.
# Usage: bash validate_bipolaris_model.sh TRAIN_DIR
set -uo pipefail

TRAIN_DIR="${1:?Usage: $0 TRAIN_DIR}"
MODEL="${AUGUSTUS_SPECIES:-bipolaris_nd90pr}"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: Activate the Conda environment first." >&2
  exit 1
fi
export AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH:-$CONDA_PREFIX/config}"

TEST="$TRAIN_DIR/ND90Pr.genes.gb.test"
[[ -s "$TEST" ]] || { echo "ERROR: Test set not found: $TEST" >&2; exit 1; }
[[ -d "$AUGUSTUS_CONFIG_PATH/species/$MODEL" ]] || {
  echo "ERROR: AUGUSTUS model not found: $MODEL" >&2
  exit 1
}

OUT="$TRAIN_DIR/ND90Pr.validation.out"
ERR="$TRAIN_DIR/ND90Pr.validation.log"

augustus --species="$MODEL" "$TEST" > "$OUT" 2> "$ERR"
status=$?
echo "AUGUSTUS exit status: $status"
echo "Prediction output: $OUT"
echo "Stderr log: $ERR"
echo "Gene records (rough count): $(grep -c '^# start gene' "$OUT" || true)"
tail -20 "$ERR" || true
exit "$status"
