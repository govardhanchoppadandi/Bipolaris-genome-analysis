#!/usr/bin/env bash
set -Eeuo pipefail

BASE="/mnt/d/Bipolaris_genomes/Bipolaris"
RNA="$BASE/NCBI_RNAseq"

mkdir -p "$RNA/sra" "$RNA/fastq" "$RNA/tmp" "$RNA/logs"

RUNS=(
    SRR15037318
    SRR15037319
    SRR15037320
    SRR15037321
)

THREADS=8

echo "========================================"
echo "BIPOLARIS RNA-SEQ DOWNLOAD AND CONVERSION"
echo "========================================"

df -h "$RNA"

for RUN in "${RUNS[@]}"; do

    echo
    echo "========================================"
    echo "PROCESSING: $RUN"
    echo "========================================"

    SRA="$RNA/sra/$RUN/$RUN.sra"
    FWD="$RNA/fastq/${RUN}_1.fastq"
    REV="$RNA/fastq/${RUN}_2.fastq"

    # 1. Download SRA if not already present
    if [[ ! -s "$SRA" ]]; then
        echo "Downloading $RUN..."
        prefetch --output-directory "$RNA/sra" "$RUN" \
            2>&1 | tee "$RNA/logs/${RUN}_prefetch.log"
    else
        echo "SRA file already exists. Skipping download."
    fi

    # 2. Validate SRA
    echo "Validating $RUN..."
    vdb-validate "$SRA" \
        2>&1 | tee "$RNA/logs/${RUN}_validate.log"

    # 3. Convert to paired FASTQ if needed
    if [[ ! -s "${FWD}.gz" || ! -s "${REV}.gz" ]]; then

        if [[ -s "$FWD" && -s "$REV" ]]; then
            echo "Uncompressed FASTQ files already exist."
        else
            echo "Converting $RUN to paired FASTQ..."

            fasterq-dump \
                --split-files \
                --threads "$THREADS" \
                --temp "$RNA/tmp" \
                --outdir "$RNA/fastq" \
                "$SRA" \
                2>&1 | tee "$RNA/logs/${RUN}_fasterq.log"
        fi

        # 4. Compress FASTQ files
        for FILE in "$FWD" "$REV"; do
            if [[ -s "$FILE" ]]; then
                echo "Compressing $(basename "$FILE")..."
                gzip -1 "$FILE"
            fi
        done
    else
        echo "Compressed FASTQ files already exist. Skipping conversion."
    fi

    # 5. Verify paired FASTQ files
    if [[ ! -s "${FWD}.gz" || ! -s "${REV}.gz" ]]; then
        echo "ERROR: Missing FASTQ pair for $RUN"
        exit 1
    fi

    echo "Checking compressed FASTQ files..."
    gzip -t "${FWD}.gz"
    gzip -t "${REV}.gz"

    echo "SUCCESS: $RUN completed."

done

echo
echo "========================================"
echo "ALL FOUR RUNS COMPLETED"
echo "========================================"

echo "FASTQ FILES:"
ls -lh "$RNA/fastq/"

echo
echo "SRA FILES:"
du -sh "$RNA/sra/"*

echo
echo "FINISHED"
