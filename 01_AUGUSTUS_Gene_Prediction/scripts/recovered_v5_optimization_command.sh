#!/usr/bin/env bash
set -Eeuo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"
SPECIES="bipolaris_multi_v5_clean"
WORK="$HOME/augustus_work/bipolaris_multi_v5_balanced_dedup"
TRAIN="$WORK/V5_raw_deduplicated_no_test_overlap.gb"
[[ -s "$TRAIN" ]] || { echo "Missing training file: $TRAIN" >&2; exit 1; }
optimize_augustus.pl --species="$SPECIES" --rounds=2 --cpus=8 "$TRAIN" 2>&1 | tee "$WORK/optimize_v5_clean.log"
