
#!/usr/bin/env bash
set -uo pipefail

# ============================================================
# BUSCO PROTEIN ASSESSMENT FOR ARI GENOMES AND BIPOLARIS
# Lineage: pleosporales_odb12
# Mode: proteins
# CPUs: 8
# ============================================================

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest

LINEAGE="pleosporales_odb12"
DOWNLOAD="$HOME/busco_downloads"
CPUS=8

# Check BUSCO
if ! command -v busco >/dev/null 2>&1; then
    echo "ERROR: BUSCO is not available in busco_latest."
    exit 1
fi

# Check lineage
if [[ ! -d "$DOWNLOAD/lineages/$LINEAGE" ]]; then
    echo "ERROR: BUSCO lineage not found:"
    echo "$DOWNLOAD/lineages/$LINEAGE"
    exit 1
fi

# ------------------------------------------------------------
# Function to process one protein collection
# ------------------------------------------------------------

run_collection() {
    local INPUT="$1"
    local OUTPUT="$2"
    local LABEL="$3"

    local REPORT="$OUTPUT/BUSCO_${LABEL}_summary.tsv"
    local MASTER="$OUTPUT/BUSCO_${LABEL}_master.log"

    mkdir -p "$OUTPUT"

    printf "Protein_file\tBUSCO_output\tStatus\n" > "$REPORT"

    echo "============================================"
    echo "Starting BUSCO: $LABEL"
    echo "Input: $INPUT"
    echo "Output: $OUTPUT"
    echo "============================================"

    if [[ ! -d "$INPUT" ]]; then
        echo "ERROR: Input folder not found: $INPUT"
        return 1
    fi

    local TOTAL=0
    local SUCCESS=0
    local SKIPPED=0
    local FAILED=0

    # Find protein FASTA files recursively.
    while IFS= read -r -d '' FILE; do

        TOTAL=$((TOTAL + 1))

        # Make a unique output name using the relative input path.
        REL="${FILE#"$INPUT"/}"
        ID="$(printf '%s' "$REL" | sed 's/[^A-Za-z0-9_.-]/_/g')"
        ID="${ID%.*}"
        OUTNAME="${LABEL}_${ID}"

        echo
        echo "Processing: $FILE" | tee -a "$MASTER"

        # Do not overwrite an existing BUSCO result.
        if [[ -d "$OUTPUT/$OUTNAME" ]]; then
            echo "SKIPPED: Result folder already exists: $OUTNAME" \
                | tee -a "$MASTER"
            printf "%s\t%s\tSKIPPED_EXISTING\n" \
                "$FILE" "$OUTNAME" >> "$REPORT"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        if busco \
            -i "$FILE" \
            -m proteins \
            -l "$LINEAGE" \
            -o "$OUTNAME" \
            --out_path "$OUTPUT" \
            -c "$CPUS" \
            --download_path "$DOWNLOAD" \
            >> "$MASTER" 2>&1; then

            echo "SUCCESS: $OUTNAME" | tee -a "$MASTER"
            printf "%s\t%s\tSUCCESS\n" \
                "$FILE" "$OUTNAME" >> "$REPORT"
            SUCCESS=$((SUCCESS + 1))

        else
            echo "FAILED: $OUTNAME" | tee -a "$MASTER"
            printf "%s\t%s\tFAILED\n" \
                "$FILE" "$OUTNAME" >> "$REPORT"
            FAILED=$((FAILED + 1))
        fi

    done < <(
        find "$INPUT" -type f \
            \( -iname "*.faa" \
            -o -iname "*.fasta" \
            -o -iname "*.fa" \
            -o -iname "*.pep" \
            -o -iname "*.fas" \) \
            -print0
    )

    echo
    echo "============================================" | tee -a "$MASTER"
    echo "FINISHED: $LABEL" | tee -a "$MASTER"
    echo "Files found: $TOTAL" | tee -a "$MASTER"
    echo "Successful: $SUCCESS" | tee -a "$MASTER"
    echo "Skipped: $SKIPPED" | tee -a "$MASTER"
    echo "Failed: $FAILED" | tee -a "$MASTER"
    echo "Report: $REPORT" | tee -a "$MASTER"
    echo "============================================" | tee -a "$MASTER"
}

# ------------------------------------------------------------
# 1. ARI GENOMES
# ------------------------------------------------------------

run_collection \
"/mnt/d/Bipolaris_genomes/ARI genomes/busco/protrin input" \
"/mnt/d/Bipolaris_genomes/ARI genomes/busco/BUSCO_RESULTS/PROTEIN" \
"ARI"

# ------------------------------------------------------------
# 2. BIPOLARIS
# ------------------------------------------------------------

run_collection \
"/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO/Protein input" \
"/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO/BUSCO_RESULTS/PROTEIN" \
"BIPOLARIS"

echo
echo "BOTH PROTEIN BUSCO COLLECTIONS FINISHED."
