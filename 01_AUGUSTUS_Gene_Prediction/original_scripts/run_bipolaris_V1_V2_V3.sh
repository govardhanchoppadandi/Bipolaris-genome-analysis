#!/usr/bin/env bash
set -Eeuo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"

GENOME="/mnt/d/Bipolaris_genomes/Bipolaris/GCA_000527765.2/GCA_000527765.2_ANU_Cvic_canu1.6_pilon5_genomic.fna"

LINEAGE="pleosporales_odb12"
CPUS=8

OUT="/mnt/d/Bipolaris_genomes/Bipolaris/FINAL_V1_V2_V3_COMPARISON"

mkdir -p "$OUT"

if [[ ! -s "$GENOME" ]]; then
    echo "ERROR: Target genome not found."
    exit 1
fi

# ------------------------------------------------
# STEP 1: AUGUSTUS GENE PREDICTION
# ------------------------------------------------

run_augustus () {
    ENV="$1"
    SPECIES="$2"
    VERSION="$3"

    DIR="$OUT/$VERSION"
    mkdir -p "$DIR"

    echo "========================================"
    echo "PREDICTING GENES: $VERSION"
    echo "Model: $SPECIES"
    echo "========================================"

    conda activate "$ENV"
    export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

    augustus \
        --species="$SPECIES" \
        --gff3=on \
        --UTR=off \
        "$GENOME" \
        > "$DIR/genes.gff3" \
        2> "$DIR/augustus.log"

    if [[ ! -s "$DIR/genes.gff3" ]]; then
        echo "ERROR: $VERSION prediction failed."
        exit 1
    fi

    echo "Predicted genes:"
    awk '$3=="gene"{n++} END{print n+0}' "$DIR/genes.gff3"
}

run_augustus bipolaris bipolaris_nd90pr V1
run_augustus busco_latest bipolaris_nd90pr_v2 V2
run_augustus busco_latest bipolaris_multi_v3_all21 V3

# ------------------------------------------------
# STEP 2: EXTRACT PROTEINS AND CDS
# ------------------------------------------------

conda activate busco_latest
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

for VERSION in V1 V2 V3; do

    DIR="$OUT/$VERSION"

    echo "Extracting proteins for $VERSION..."

    gffread "$DIR/genes.gff3" \
        -g "$GENOME" \
        -y "$DIR/proteins.faa" \
        -x "$DIR/cds.fna"

    if [[ ! -s "$DIR/proteins.faa" ]]; then
        echo "ERROR: $VERSION protein FASTA is empty."
        exit 1
    fi

    echo "Proteins: $(grep -c '^>' "$DIR/proteins.faa")"
    echo "CDS: $(grep -c '^>' "$DIR/cds.fna")"
done

# ------------------------------------------------
# STEP 3: BUSCO ON PREDICTED PROTEINS
# ------------------------------------------------

for VERSION in V1 V2 V3; do

    DIR="$OUT/$VERSION"

    echo "========================================"
    echo "BUSCO PROTEIN ANALYSIS: $VERSION"
    echo "========================================"

    busco \
        -i "$DIR/proteins.faa" \
        -m proteins \
        -l "$LINEAGE" \
        -o "${VERSION}_BUSCO" \
        --out_path "$DIR" \
        --cpu "$CPUS"

done

# ------------------------------------------------
# STEP 4: FINAL COMPARISON
# ------------------------------------------------

SUMMARY="$OUT/FINAL_COMPARISON.txt"

{
    echo "BIPOLARIS AUGUSTUS V1 vs V2 vs V3"
    echo "Target: GCA_000527765.2"
    echo "BUSCO lineage: $LINEAGE"
    echo "BUSCO mode: proteins"
    echo

    for VERSION in V1 V2 V3; do

        DIR="$OUT/$VERSION"

        echo "----------------------------------------"
        echo "$VERSION"
        echo "----------------------------------------"

        echo "Genes:"
        awk '$3=="gene"{n++} END{print n+0}' "$DIR/genes.gff3"

        echo "Proteins:"
        grep -c '^>' "$DIR/proteins.faa"

        echo "CDS:"
        grep -c '^>' "$DIR/cds.fna"

        REPORT=$(find "$DIR" -name "short_summary*.txt" -type f | head -n 1)

        if [[ -n "$REPORT" ]]; then
            cat "$REPORT"
        else
            echo "BUSCO summary not found."
        fi

        echo
    done
} > "$SUMMARY"

echo
echo "========================================"
echo "ALL PREDICTIONS AND BUSCO RUNS COMPLETED"
echo "========================================"
echo "Results: $OUT"
echo "Comparison: $SUMMARY"
cat "$SUMMARY"

