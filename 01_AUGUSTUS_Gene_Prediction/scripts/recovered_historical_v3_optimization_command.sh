#!/usr/bin/env bash
set -Eeuo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"
SPECIES="bipolaris_multi_v3"
WORK="$HOME/augustus_work/bipolaris_multi_v3_verified_20261102_093515"
TRAIN="$WORK/bipolaris_multi_v3.all_references.gb.train"
METAPARS="$AUGUSTUS_CONFIG_PATH/species/$SPECIES/${SPECIES}_metapars.cfg"
[[ -s "$TRAIN" && -s "$METAPARS" ]] || { echo "Historical V3 training or metaparameter file missing." >&2; exit 1; }
optimize_augustus.pl --species="$SPECIES" --metapars="$METAPARS" --rounds=2 --cpus=8 --cleanup=1 "$TRAIN"
