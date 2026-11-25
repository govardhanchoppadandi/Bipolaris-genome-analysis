
#!/usr/bin/env bash

set -uo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate busco_latest

export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

ROOT="/mnt/d/Bipolaris_genomes/Bipolaris"
MODEL="bipolaris_nd90pr_v2"

MODEL_DIR="$AUGUSTUS_CONFIG_PATH/species/$MODEL"

# Check programs and model
for program in augustus gffread; do
    command -v "$program" >/dev/null 2>&1 || {
        echo "ERROR: $program not found."
        exit 1
    }
done

[ -d "$ROOT" ] || {
    echo "ERROR: Genome directory not found: $ROOT"
    exit 1
}

[ -d "$MODEL_DIR" ] || {
    echo "ERROR: Trained AUGUSTUS model not found: $MODEL_DIR"
    exit 1
}

PROCESSED=0
SKIPPED_ANNOTATED=0
SKIPPED_NO_GENOME=0
SKIPPED_EXISTING=0
FAILED=0

echo "=========================================="
echo "BIPOLARIS AUGUSTUS ANNOTATION PIPELINE"
echo "=========================================="
echo "AUGUSTUS model: $MODEL"
echo "Model directory: $MODEL_DIR"
echo "Input directory: $ROOT"
echo "AUGUSTUS executable: $(command -v augustus)"
echo "gffread executable: $(command -v gffread)"
echo "=========================================="

for DIR in "$ROOT"/*/; do

    [ -d "$DIR" ] || continue

    ACCESSION=$(basename "$DIR")

    echo
    echo "------------------------------------------"
    echo "Genome accession: $ACCESSION"
    echo "MODEL USED: $MODEL"
    echo "------------------------------------------"

    # Skip already annotated genomes
    if [ -f "$DIR/protein.faa" ]; then
        echo "SKIPPED: Existing protein.faa found."
        ((SKIPPED_ANNOTATED+=1))
        continue
    fi

    # Find genome FASTA, prioritizing genomic.fna
    GENOME=""

    for PATTERN in \
        '*_genomic.fna' \
        '*.fna' \
        '*.fasta' \
        '*.fa' \
        '*.fas'
    do
        for CANDIDATE in "$DIR"/$PATTERN; do

            [ -f "$CANDIDATE" ] || continue

            BASE=$(basename "$CANDIDATE")
            LOWER=${BASE,,}

            # Exclude non-genome sequence files
            case "$LOWER" in
                *protein*|*cds*|*rna*|*transcript*|*mrna*|*pep*)
                    continue
                    ;;
            esac

            # Require FASTA headers
            if ! grep -q '^>' "$CANDIDATE"; then
                continue
            fi

            GENOME="$CANDIDATE"
            break 2
        done
    done

    if [ -z "$GENOME" ]; then
        echo "SKIPPED: No suitable genome FASTA found."
        ((SKIPPED_NO_GENOME+=1))
        continue
    fi

    GFF="$DIR/augustus_genes.gff3"
    PROTEIN="$DIR/augustus_proteins.faa"
    CDS="$DIR/augustus_cds.fna"
    LOG="$DIR/augustus_prediction.log"

    # Do not overwrite previous results
    if [ -e "$GFF" ] || [ -e "$PROTEIN" ] || [ -e "$CDS" ]; then
        echo "SKIPPED: Previous AUGUSTUS output exists."
        ((SKIPPED_EXISTING+=1))
        continue
    fi

    GFF_TMP="$DIR/.augustus_genes.tmp.$$"
    PROTEIN_TMP="$DIR/.augustus_proteins.tmp.$$"
    CDS_TMP="$DIR/.augustus_cds.tmp.$$"

    echo "Genome FASTA: $GENOME"
    echo "MODEL USED: $MODEL"
    echo "Starting gene prediction..."

    if ! augustus \
        --species="$MODEL" \
        --gff3=on \
        "$GENOME" \
        > "$GFF_TMP" 2> "$LOG"; then

        echo "FAILED: AUGUSTUS prediction."
        rm -f "$GFF_TMP" "$PROTEIN_TMP" "$CDS_TMP"
        ((FAILED+=1))
        continue
    fi

    echo "AUGUSTUS prediction completed."
    echo "Extracting protein and CDS sequences..."

    if ! gffread "$GFF_TMP" \
        -g "$GENOME" \
        -y "$PROTEIN_TMP" \
        -x "$CDS_TMP" >> "$LOG" 2>&1; then

        echo "FAILED: gffread extraction."
        rm -f "$GFF_TMP" "$PROTEIN_TMP" "$CDS_TMP"
        ((FAILED+=1))
        continue
    fi

    GENE_COUNT=$(awk -F '\t' '$3=="gene"{n++} END{print n+0}' "$GFF_TMP")
    PROTEIN_COUNT=$(grep -c '^>' "$PROTEIN_TMP" || true)
    CDS_COUNT=$(grep -c '^>' "$CDS_TMP" || true)

    if [ "$GENE_COUNT" -eq 0 ] || \
       [ "$PROTEIN_COUNT" -eq 0 ] || \
       [ "$CDS_COUNT" -eq 0 ]; then

        echo "FAILED: Empty gene, protein or CDS output."
        rm -f "$GFF_TMP" "$PROTEIN_TMP" "$CDS_TMP"
        ((FAILED+=1))
        continue
    fi

    mv "$GFF_TMP" "$GFF"
    mv "$PROTEIN_TMP" "$PROTEIN"
    mv "$CDS_TMP" "$CDS"

    echo "SUCCESS: $ACCESSION"
    echo "Model: $MODEL"
    echo "Genes: $GENE_COUNT"
    echo "Proteins: $PROTEIN_COUNT"
    echo "CDS: $CDS_COUNT"
    echo "Saved inside: $DIR"

    ((PROCESSED+=1))

done

echo
echo "=========================================="
echo "FINAL ANNOTATION SUMMARY"
echo "=========================================="
echo "Model used: $MODEL"
echo "Successfully processed: $PROCESSED"
echo "Skipped (protein.faa exists): $SKIPPED_ANNOTATED"
echo "Skipped (no genome FASTA): $SKIPPED_NO_GENOME"
echo "Skipped (previous AUGUSTUS outputs): $SKIPPED_EXISTING"
echo "Failed: $FAILED"
echo "=========================================="
