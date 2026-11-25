
#!/usr/bin/env bash
set -Eeuo pipefail

# ==========================================================
# ARI GENOMES - AUGUSTUS V2 PREDICTION
# ==========================================================

ROOT="/mnt/d/Bipolaris_genomes/ARI genomes"
MODEL="bipolaris_nd90pr_v2"

LOGDIR="$ROOT/V2_AUGUSTUS_LOGS"
REPORT="$ROOT/V2_AUGUSTUS_prediction_report.tsv"

mkdir -p "$LOGDIR"

# Activate Conda
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest

# Configure AUGUSTUS
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

# Check required software and model
for CMD in augustus gffread; do
    if ! command -v "$CMD" >/dev/null 2>&1; then
        echo "ERROR: $CMD is not installed or not in PATH."
        exit 1
    fi
done

MODEL_DIR="$AUGUSTUS_CONFIG_PATH/species/$MODEL"

if [[ ! -d "$MODEL_DIR" ]]; then
    echo "ERROR: V2 model not found:"
    echo "$MODEL_DIR"
    exit 1
fi

if [[ ! -d "$ROOT" ]]; then
    echo "ERROR: ARI genome directory not found: $ROOT"
    exit 1
fi

# Initialize report
printf "Accession\tGenome\tStatus\tGenes\tProteins\tCDS\tLog\n" > "$REPORT"

MASTER="$LOGDIR/V2_master.log"

echo "==============================================" | tee "$MASTER"
echo "ARI GENOMES AUGUSTUS V2 PREDICTION" | tee -a "$MASTER"
echo "Date: $(date)" | tee -a "$MASTER"
echo "Model: $MODEL" | tee -a "$MASTER"
echo "Root: $ROOT" | tee -a "$MASTER"
echo "==============================================" | tee -a "$MASTER"

TOTAL=0
SUCCESS=0
SKIPPED=0
FAILED=0

# Process immediate accession folders only.
while IFS= read -r -d '' DIR; do

    ACCESSION="$(basename "$DIR")"

    # Explicitly ignore BUSCO directory
    if [[ "${ACCESSION,,}" == "busco" ]]; then
        echo "SKIP: BUSCO folder" | tee -a "$MASTER"
        continue
    fi

    # Ignore our own output folders
    if [[ "$ACCESSION" == "V2_AUGUSTUS_LOGS" ]]; then
        continue
    fi

    TOTAL=$((TOTAL + 1))

    LOG="$LOGDIR/${ACCESSION}.log"

    echo | tee -a "$MASTER"
    echo "Processing: $ACCESSION" | tee -a "$MASTER"

    # Skip if existing protein FASTA is present
    if [[ -f "$DIR/protein.faa" || -f "$DIR/proteins.faa" ]]; then
        echo "SKIPPED: Existing protein FASTA" | tee -a "$MASTER"
        printf "%s\tNA\tSKIPPED_EXISTING_PROTEIN\t0\t0\t0\t%s\n" \
            "$ACCESSION" "$LOG" >> "$REPORT"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Avoid overwriting partial or previous prediction outputs
    if [[ -e "$DIR/genes.gff3" || -e "$DIR/cds.fna" ]]; then
        echo "SKIPPED: Existing genes.gff3 or cds.fna" | tee -a "$MASTER"
        printf "%s\tNA\tSKIPPED_EXISTING_OUTPUT\t0\t0\t0\t%s\n" \
            "$ACCESSION" "$LOG" >> "$REPORT"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Locate genome FASTA files directly inside accession folder
    mapfile -d '' GENOMES < <(
        find "$DIR" -maxdepth 1 -type f \
            \( -iname "*.fna" -o -iname "*.fasta" -o -iname "*.fa" \) \
            ! -iname "*protein*" \
            ! -iname "*proteins*" \
            ! -iname "*cds*" \
            ! -iname "*transcript*" \
            ! -iname "*rna*" \
            -print0
    )

    if [[ ${#GENOMES[@]} -eq 0 ]]; then
        echo "SKIPPED: No genome FASTA found" | tee -a "$MASTER"
        printf "%s\tNA\tSKIPPED_NO_GENOME\t0\t0\t0\t%s\n" \
            "$ACCESSION" "$LOG" >> "$REPORT"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # If multiple candidates, prefer a unique genomic filename
    if [[ ${#GENOMES[@]} -gt 1 ]]; then
        GENOMIC=()

        for F in "${GENOMES[@]}"; do
            if [[ "$(basename "$F")" =~ genomic ]]; then
                GENOMIC+=("$F")
            fi
        done

        if [[ ${#GENOMIC[@]} -eq 1 ]]; then
            GENOME="${GENOMIC[0]}"
        else
            echo "SKIPPED: Multiple ambiguous genome FASTA files" | tee -a "$MASTER"
            printf "%s\tNA\tSKIPPED_AMBIGUOUS_GENOME\t0\t0\t0\t%s\n" \
                "$ACCESSION" "$LOG" >> "$REPORT"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi
    else
        GENOME="${GENOMES[0]}"
    fi

    echo "Genome: $GENOME" | tee -a "$MASTER"

    TMPGFF="$DIR/.genes_V2.tmp.gff3"
    TMPPROT="$DIR/.proteins_V2.tmp.faa"
    TMPCDS="$DIR/.cds_V2.tmp.fna"

    rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"

    # Run AUGUSTUS
    if ! augustus \
        --species="$MODEL" \
        --gff3=on \
        --UTR=off \
        "$GENOME" \
        > "$TMPGFF" 2> "$LOG"; then

        echo "FAILED: AUGUSTUS prediction" | tee -a "$MASTER"
        printf "%s\t%s\tFAILED_AUGUSTUS\t0\t0\t0\t%s\n" \
            "$ACCESSION" "$GENOME" "$LOG" >> "$REPORT"

        rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Extract protein and CDS sequences
    if ! gffread "$TMPGFF" \
        -g "$GENOME" \
        -y "$TMPPROT" \
        -x "$TMPCDS" \
        2>> "$LOG"; then

        echo "FAILED: gffread extraction" | tee -a "$MASTER"
        printf "%s\t%s\tFAILED_GFFREAD\t0\t0\t0\t%s\n" \
            "$ACCESSION" "$GENOME" "$LOG" >> "$REPORT"

        rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Count output sequences
    GENES=$(grep -c '^>' "$TMPPROT" || true)
    PROTEINS=$(grep -c '^>' "$TMPPROT" || true)
    CDS=$(grep -c '^>' "$TMPCDS" || true)

    # Check generated files
    if [[ ! -s "$TMPGFF" || ! -s "$TMPPROT" || ! -s "$TMPCDS" ]]; then
        echo "FAILED: Empty output file" | tee -a "$MASTER"
        printf "%s\t%s\tFAILED_EMPTY_OUTPUT\t%s\t%s\t%s\t%s\n" \
            "$ACCESSION" "$GENOME" "$GENES" "$PROTEINS" "$CDS" "$LOG" \
            >> "$REPORT"

        rm -f "$TMPGFF" "$TMPPROT" "$TMPCDS"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Save final outputs inside accession folder
    mv "$TMPGFF" "$DIR/genes.gff3"
    mv "$TMPPROT" "$DIR/proteins.faa"
    mv "$TMPCDS" "$DIR/cds.fna"

    echo "SUCCESS: $ACCESSION | Genes: $GENES | Proteins: $PROTEINS | CDS: $CDS" \
        | tee -a "$MASTER"

    printf "%s\t%s\tSUCCESS\t%s\t%s\t%s\t%s\n" \
        "$ACCESSION" "$GENOME" "$GENES" "$PROTEINS" "$CDS" "$LOG" \
        >> "$REPORT"

    SUCCESS=$((SUCCESS + 1))

done < <(
    find "$ROOT" -mindepth 1 -maxdepth 1 -type d \
        ! -iname "busco" -print0
)

echo | tee -a "$MASTER"
echo "==============================================" | tee -a "$MASTER"
echo "ARI V2 PREDICTION FINISHED" | tee -a "$MASTER"
echo "Folders examined: $TOTAL" | tee -a "$MASTER"
echo "Successful: $SUCCESS" | tee -a "$MASTER"
echo "Skipped: $SKIPPED" | tee -a "$MASTER"
echo "Failed: $FAILED" | tee -a "$MASTER"
echo "Report: $REPORT" | tee -a "$MASTER"
echo "==============================================" | tee -a "$MASTER"
