
#!/usr/bin/env bash
set -uo pipefail
export LC_ALL=C

# ============================================================
# BIPOLARIS ALL-GENOME PREDICTION USING AUGUSTUS V2
#
# Model: bipolaris_nd90pr_v2
# Input: all accession folders directly inside Bipolaris/
#
# Skip:
#   - protein.faa already exists
#   - proteins.faa already exists
#   - no whole-genome FASTA
#   - ambiguous multiple genome FASTAs
#
# Output in each accession folder:
#   genes.gff3
#   proteins.faa
#   cds.fna
#
# ============================================================

# ---------------- ENVIRONMENT ----------------

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest

SPECIES="bipolaris_nd90pr_v2"
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

# ---------------- INPUT ROOT ----------------

ROOT="/mnt/d/Bipolaris_genomes/Bipolaris"

# ---------------- LOGS AND REPORT ----------------

REPORT="$ROOT/V2_gene_prediction_report.tsv"
LOGDIR="$ROOT/V2_AUGUSTUS_LOGS"
MASTERLOG="$LOGDIR/V2_master.log"

mkdir -p "$LOGDIR"

# ---------------- CHECK SOFTWARE ----------------

for program in augustus gffread find; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program"
        exit 1
    fi
done

if [[ ! -d "$AUGUSTUS_CONFIG_PATH/species/$SPECIES" ]]; then
    echo "ERROR: V2 model not found:"
    echo "$AUGUSTUS_CONFIG_PATH/species/$SPECIES"
    exit 1
fi

if [[ ! -d "$ROOT" ]]; then
    echo "ERROR: Input directory not found: $ROOT"
    exit 1
fi

# ---------------- INITIALIZE REPORT ----------------

printf 'Genome_Folder\tGenome_FASTA\tModel\tGenes\tProteins\tCDS\tStatus\tReason\n' \
    > "$REPORT"

: > "$MASTERLOG"

# ---------------- COUNTERS ----------------

TOTAL=0
COMPLETED=0
SKIPPED=0
FAILED=0
NO_FASTA=0
AMBIGUOUS=0
CONFLICT=0

echo
echo "================================================"
echo "BIPOLARIS V2 BATCH GENE PREDICTION"
echo "================================================"
echo "Model: $SPECIES"
echo "Input: $ROOT"
echo "Output: Individual accession folders"
echo "================================================"

# ---------------- PROCESS ACCESSION FOLDERS ----------------

while IFS= read -r -d '' DIR; do

    TOTAL=$((TOTAL+1))

    NAME="$(basename "$DIR")"

    echo
    echo "------------------------------------------------"
    echo "ACCESSION: $NAME"
    echo "------------------------------------------------"

    # Skip if existing protein FASTA is present.
    if [[ -e "$DIR/protein.faa" || -e "$DIR/proteins.faa" ]]; then

        echo "SKIPPED: Existing protein FASTA found."

        printf '%s\t\t%s\t\t\t\tSKIPPED\tExisting protein.faa or proteins.faa\n' \
            "$NAME" "$SPECIES" >> "$REPORT"

        SKIPPED=$((SKIPPED+1))
        continue
    fi

    # Never overwrite existing prediction outputs.
    if [[ -e "$DIR/genes.gff3" || -e "$DIR/cds.fna" ]]; then

        echo "SKIPPED: Existing prediction output conflict."

        printf '%s\t\t%s\t\t\t\tOUTPUT_CONFLICT\tExisting genes.gff3 or cds.fna\n' \
            "$NAME" "$SPECIES" >> "$REPORT"

        CONFLICT=$((CONFLICT+1))
        continue
    fi

    # ---------------- FIND GENOME FASTA ----------------

    shopt -s nullglob nocaseglob

    ALL_FILES=(
        "$DIR"/*.fna
        "$DIR"/*.fasta
        "$DIR"/*.fa
    )

    shopt -u nullglob nocaseglob

    CANDIDATES=()

    for FILE in "${ALL_FILES[@]}"; do

        [[ -f "$FILE" ]] || continue

        FILE_NAME="$(basename "$FILE")"

        # Exclude protein, CDS, RNA and transcript sequences.
        if [[ "$FILE_NAME" =~ (protein|proteins|cds|transcript|rna) ]]; then
            continue
        fi

        CANDIDATES+=("$FILE")

    done

    # No usable genome FASTA.
    if (( ${#CANDIDATES[@]} == 0 )); then

        echo "NO FASTA: No whole-genome FASTA found."

        printf '%s\t\t%s\t\t\t\tNO_FASTA\tNo suitable genome FASTA\n' \
            "$NAME" "$SPECIES" >> "$REPORT"

        NO_FASTA=$((NO_FASTA+1))
        continue
    fi

    # If several FASTAs exist, use a unique genomic file if available.
    if (( ${#CANDIDATES[@]} > 1 )); then

        GENOMIC=()

        for FILE in "${CANDIDATES[@]}"; do
            if [[ "$(basename "$FILE")" =~ genomic ]]; then
                GENOMIC+=("$FILE")
            fi
        done

        if (( ${#GENOMIC[@]} == 1 )); then
            CANDIDATES=("${GENOMIC[@]}")
        else

            echo "AMBIGUOUS: Multiple genome FASTA files found."
            printf '%s\n' "${CANDIDATES[@]}"

            printf '%s\t\t%s\t\t\t\tAMBIGUOUS_FASTA\tMultiple possible genome FASTAs\n' \
                "$NAME" "$SPECIES" >> "$REPORT"

            AMBIGUOUS=$((AMBIGUOUS+1))
            continue
        fi
    fi

    GENOME="${CANDIDATES[0]}"

    if [[ ! -s "$GENOME" ]]; then
        echo "FAILED: Genome FASTA is empty."
        FAILED=$((FAILED+1))
        continue
    fi

    echo "Genome FASTA: $GENOME"

    # ---------------- OUTPUT FILES ----------------

    GFF="$DIR/genes.gff3"
    PROTEIN="$DIR/proteins.faa"
    CDS="$DIR/cds.fna"

    TMPGFF="$DIR/genes.gff3.tmp.$$"
    TMPPROT="$DIR/proteins.faa.tmp.$$"
    TMPCDS="$DIR/cds.fna.tmp.$$"

    LOG="$LOGDIR/${NAME}_V2.log"

    # ---------------- AUGUSTUS PREDICTION ----------------

    echo "Running AUGUSTUS V2..."

    if augustus \
        --species="$SPECIES" \
        --gff3=on \
        --UTR=off \
        "$GENOME" \
        > "$TMPGFF" 2> "$LOG"; then

        echo "AUGUSTUS completed."

    else

        echo "FAILED: AUGUSTUS prediction."

        printf '%s\t%s\t%s\t\t\t\tFAILED\tAUGUSTUS error; see log\n' \
            "$NAME" "$GENOME" "$SPECIES" >> "$REPORT"

        rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"

        echo "AUGUSTUS error for $NAME" >> "$MASTERLOG"

        FAILED=$((FAILED+1))
        continue
    fi

    if [[ ! -s "$TMPGFF" ]]; then

        echo "FAILED: Empty GFF3 output."

        rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"

        FAILED=$((FAILED+1))
        continue
    fi

    # ---------------- PROTEIN AND CDS EXTRACTION ----------------

    echo "Extracting proteins and CDS using gffread..."

    if gffread "$TMPGFF" \
        -g "$GENOME" \
        -y "$TMPPROT" \
        -x "$TMPCDS" \
        2>> "$LOG" &&
        [[ -s "$TMPPROT" && -s "$TMPCDS" ]]; then

        echo "Sequence extraction completed."

    else

        echo "FAILED: gffread extraction."

        printf '%s\t%s\t%s\t\t\t\tFAILED\tgffread extraction error; see log\n' \
            "$NAME" "$GENOME" "$SPECIES" >> "$REPORT"

        rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"

        FAILED=$((FAILED+1))
        continue
    fi

    # ---------------- COUNT RESULTS ----------------

    GENES=$(awk -F '\t' '$3=="gene"{n++} END{print n+0}' "$TMPGFF")

    PROTEINS=$(grep -c '^>' "$TMPPROT" || true)
    CDS_COUNT=$(grep -c '^>' "$TMPCDS" || true)

    if (( PROTEINS == 0 || CDS_COUNT == 0 )); then

        echo "FAILED: No usable proteins or CDS."

        rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"

        FAILED=$((FAILED+1))
        continue
    fi

    # ---------------- SAVE RESULTS ----------------

    mv "$TMPGFF" "$GFF"
    mv "$TMPPROT" "$PROTEIN"
    mv "$TMPCDS" "$CDS"

    echo
    echo "SUCCESS: $NAME"
    echo "Genes: $GENES"
    echo "Proteins: $PROTEINS"
    echo "CDS: $CDS_COUNT"

    printf '%s\t%s\t%s\t%s\t%s\t%s\tCOMPLETED\tPrediction and extraction successful\n' \
        "$NAME" "$GENOME" "$SPECIES" \
        "$GENES" "$PROTEINS" "$CDS_COUNT" >> "$REPORT"

    COMPLETED=$((COMPLETED+1))

done < <(
    find "$ROOT" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        ! -name "V2_AUGUSTUS_LOGS" \
        -print0 | sort -z
)

# ---------------- FINAL SUMMARY ----------------

{
    echo
    echo "================================================"
    echo "V2 BATCH PREDICTION COMPLETED"
    echo "================================================"
    echo "Input directory: $ROOT"
    echo "AUGUSTUS model: $SPECIES"
    echo "Total folders: $TOTAL"
    echo "Completed: $COMPLETED"
    echo "Skipped existing proteins: $SKIPPED"
    echo "Failed: $FAILED"
    echo "No genome FASTA: $NO_FASTA"
    echo "Ambiguous genome FASTA: $AMBIGUOUS"
    echo "Output conflicts: $CONFLICT"
    echo
    echo "REPORT: $REPORT"
    echo "LOGS: $LOGDIR"
    echo "================================================"
} | tee -a "$MASTERLOG"

echo
echo "V2 prediction batch finished."
