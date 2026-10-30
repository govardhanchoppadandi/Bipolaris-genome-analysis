#!/usr/bin/env bash
set -Eeuo pipefail

BASE="/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v busco >/dev/null || {
    echo "Activate busco_latest first:"
    echo "source ~/miniconda3/etc/profile.d/conda.sh"
    echo "conda activate busco_latest"
    exit 1
}

echo "Bipolaris BUSCO workflow"
busco --version
echo "Executable: $(command -v busco)"

mkdir -p "$BASE/BUSCO_RESULTS/GENOME" "$BASE/BUSCO_RESULTS/PROTEIN" \
 "$BASE/BUSCO_RESULTS/FAILED" "$BASE/BUSCO_RESULTS/LOGS" \
 "$BASE/BUSCO_RESULTS/DOWNLOADS" "$BASE/results"

bash "$SCRIPT_DIR/02_download_lineage.sh"

shopt -s nullglob
GENOMES=("$BASE/Genome input"/*.fasta)
PROTEINS=("$BASE/Protein input"/*.faa)
echo "Genome FASTA files: ${#GENOMES[@]}"
echo "Protein FAA files: ${#PROTEINS[@]}"
((${#GENOMES[@]} > 0)) || { echo "No genome .fasta files"; exit 1; }
((${#PROTEINS[@]} > 0)) || { echo "No protein .faa files"; exit 1; }

bash "$SCRIPT_DIR/03_run_genomes.sh"
bash "$SCRIPT_DIR/04_run_proteins.sh"
python "$SCRIPT_DIR/05_summarize_results.py"

echo "Genome TXT summaries: $(find "$BASE/BUSCO_RESULTS/GENOME" -type f -name 'short_summary*.txt' | wc -l)"
echo "Protein TXT summaries: $(find "$BASE/BUSCO_RESULTS/PROTEIN" -type f -name 'short_summary*.txt' | wc -l)"
echo "Failed folders: $(find "$BASE/BUSCO_RESULTS/FAILED" -mindepth 1 -maxdepth 1 -type d | wc -l)"
echo "Reports:"
ls -lh "$BASE/results"
echo "Pipeline finished."
