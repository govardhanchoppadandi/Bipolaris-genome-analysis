#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

OUT="results/model_validation"
mkdir -p "$OUT"

REPORT="$OUT/prediction_output_integrity.tsv"

printf "Model\tGenes\tTranscripts\tCDS_records\tProtein_records\tStatus\n" > "$REPORT"

models=(
  bipolaris_nd90pr
  bipolaris_nd90pr_v2
  bipolaris_multi_v3_all21
  bipolaris_multi_v4_advanced
  bipolaris_multi_v5_clean
  bipolaris_multi_v5_deduplicated
)

for model in "${models[@]}"; do
    gff="test_predictions/${model}.gff3"
    cds="test_predictions/cds/${model}.fna"
    proteins="test_predictions/proteins/${model}.faa"

    for file in "$gff" "$cds" "$proteins"; do
        if [[ ! -s "$file" ]]; then
            echo "ERROR: Missing or empty file: $file" >&2
            exit 1
        fi
    done

    genes=$(awk -F '\t' '!/^#/ && $3=="gene"{n++} END{print n+0}' "$gff")
    transcripts=$(awk -F '\t' '!/^#/ && ($3=="transcript" || $3=="mRNA"){n++} END{print n+0}' "$gff")
    cds_count=$(grep -c '^>' "$cds" || true)
    protein_count=$(grep -c '^>' "$proteins" || true)

    status="PASS"
    if [[ "$genes" -ne "$transcripts" ||
          "$genes" -ne "$cds_count" ||
          "$genes" -ne "$protein_count" ]]; then
        status="FAIL"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$model" "$genes" "$transcripts" "$cds_count" \
        "$protein_count" "$status" >> "$REPORT"
done

echo "===== PREDICTION OUTPUT INTEGRITY ====="
column -t -s $'\t' "$REPORT"

if grep -q $'\tFAIL$' "$REPORT"; then
    echo "QC RESULT: FAIL"
    exit 1
fi

echo
echo "QC RESULT: PASS — all six models have matching counts."
