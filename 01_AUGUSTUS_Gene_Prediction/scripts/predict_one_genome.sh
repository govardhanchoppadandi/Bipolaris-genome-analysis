#!/usr/bin/env bash
set -Eeuo pipefail
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 AUGUSTUS_SPECIES genome.fna output_directory" >&2
  exit 2
fi
SPECIES="$1"; GENOME="$2"; OUT="$3"
command -v augustus >/dev/null
command -v gffread >/dev/null
[[ -s "$GENOME" ]] || { echo "Genome FASTA missing/empty: $GENOME" >&2; exit 1; }
mkdir -p "$OUT"
export AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH:-${CONDA_PREFIX:?Activate AUGUSTUS Conda environment}/config}"
augustus --species="$SPECIES" --gff3=on --UTR=off "$GENOME" > "$OUT/genes.gff3"
gffread "$OUT/genes.gff3" -g "$GENOME" -y "$OUT/proteins.faa" -x "$OUT/cds.fna"
test -s "$OUT/genes.gff3" && test -s "$OUT/proteins.faa" && test -s "$OUT/cds.fna"
echo "Finished: $OUT"
