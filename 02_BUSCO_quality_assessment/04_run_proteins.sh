#!/usr/bin/env bash
set -Eeuo pipefail
BASE="/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO"
LINEAGE="pleosporales_odb12"
THREADS="${THREADS:-8}"
DOWNLOAD_PATH="$HOME/busco_downloads"
INPUT="$BASE/Protein input"
OUTROOT="$BASE/BUSCO_RESULTS/PROTEIN"
FAILED="$BASE/BUSCO_RESULTS/FAILED"
LOGS="$BASE/BUSCO_RESULTS/LOGS"

mkdir -p "$OUTROOT" "$FAILED" "$LOGS"
command -v busco >/dev/null || { echo "Activate busco_latest first"; exit 1; }
test -d "$DOWNLOAD_PATH/lineages/$LINEAGE" || { echo "Lineage missing; run 02_download_lineage.sh"; exit 1; }

shopt -s nullglob
FILES=("$INPUT"/*.faa)
echo "Protein files found: ${#FILES[@]}"
((${#FILES[@]} > 0)) || { echo "No .faa files found"; exit 1; }

for FILE in "${FILES[@]}"; do
    NAME="$(basename "$FILE" .faa)"
    OUT="$OUTROOT/$NAME"
    LOG="$LOGS/PROTEIN_${NAME}.log"

    if [[ -d "$OUT" ]] &&
       find "$OUT" -maxdepth 1 -type f -name 'short_summary*.txt' -size +0c | grep -q . &&
       find "$OUT" -maxdepth 1 -type f -name 'short_summary*.json' -size +0c | grep -q .; then
        echo "SKIP completed protein: $NAME"
        continue
    fi

    if [[ -d "$OUT" ]]; then
        mv "$OUT" "$FAILED/PROTEIN_${NAME}_incomplete_$(date +%Y%m%d_%H%M%S)"
    fi

    echo "RUNNING PROTEIN: $NAME"
    if busco -i "$FILE" -m proteins -l "$LINEAGE" -o "$NAME" \
        --out_path "$OUTROOT" -c "$THREADS" --download_path "$DOWNLOAD_PATH" \
        >"$LOG" 2>&1; then

        TXT="$(find "$OUT" -type f -name 'short_summary*.txt' -print -quit || true)"
        JSON="$(find "$OUT" -type f -name 'short_summary*.json' -print -quit || true)"
        if [[ -n "$TXT" && -n "$JSON" ]]; then
            TMP="$(mktemp -d "$OUT/.summary.XXXXXX")"
            cp "$TXT" "$JSON" "$TMP/"
            find "$OUT" -mindepth 1 -maxdepth 1 ! -path "$TMP" -exec rm -rf -- {} +
            mv "$TMP"/* "$OUT/"
            rmdir "$TMP"
            echo "SUCCESS: $NAME (TXT + JSON retained)"
        else
            [[ -d "$OUT" ]] && mv "$OUT" "$FAILED/PROTEIN_${NAME}_missing_summary_$(date +%Y%m%d_%H%M%S)"
            echo "FAILED: missing summary; log: $LOG"
        fi
    else
        [[ -d "$OUT" ]] && mv "$OUT" "$FAILED/PROTEIN_${NAME}_error_$(date +%Y%m%d_%H%M%S)"
        echo "FAILED: $NAME; log: $LOG"
    fi
done
echo "Protein stage finished."
