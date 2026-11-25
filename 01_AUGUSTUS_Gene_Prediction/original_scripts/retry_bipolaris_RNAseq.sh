#!/usr/bin/env bash
set -Eeuo pipefail

conda activate sra_download

OUT="/mnt/d/Bipolaris_genomes/Bipolaris/NCBI_RNAseq_V4"

mkdir -p "$OUT/metadata" "$OUT/sra" "$OUT/fastq" "$OUT/temp"

echo "Retrieving NCBI SRA metadata in smaller requests..."

PROJECTS=(
PRJNA1131952
PRJNA1162285
PRJNA847724
PRJNA743515
PRJNA492913
PRJNA1054994
)

for PROJECT in "${PROJECTS[@]}"; do

    echo "Searching $PROJECT"

    for ATTEMPT in 1 2 3 4 5; do

        if esearch -db sra \
            -query "\"$PROJECT\"[BioProject] AND \"RNA-Seq\"[Strategy]" |
            efetch -format runinfo \
            > "$OUT/metadata/${PROJECT}.csv.tmp"; then

            mv "$OUT/metadata/${PROJECT}.csv.tmp" \
               "$OUT/metadata/${PROJECT}.csv"

            break
        fi

        echo "Retry $ATTEMPT for $PROJECT"
        sleep 10
    done

    if [[ ! -s "$OUT/metadata/${PROJECT}.csv" ]]; then
        echo "WARNING: Could not retrieve $PROJECT"
    fi

    sleep 2
done

# Combine project metadata
python - "$OUT/metadata" <<'PY'
import csv, glob, os, sys

directory = sys.argv[1]
records = {}

for filename in glob.glob(directory + "/PRJNA*.csv"):
    try:
        with open(filename, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("LibraryStrategy", "").upper() == "RNA-SEQ" and \
                   row.get("LibrarySource", "").upper() == "TRANSCRIPTOMIC":
                    run = row.get("Run", "").strip()
                    if run:
                        records[run] = row
    except Exception as e:
        print("Problem reading", filename, e)

if records:
    with open(directory + "/Bipolaris_RNAseq_manifest.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(next(iter(records.values())).keys()))
        writer.writeheader()
        writer.writerows(records.values())

    with open(directory + "/run_accessions.txt", "w") as f:
        for run in sorted(records):
            f.write(run + "\n")

print("Unique RNA-seq runs:", len(records))
PY

RUNFILE="$OUT/metadata/run_accessions.txt"

if [[ ! -s "$RUNFILE" ]]; then
    echo "No usable run list generated. Check the metadata."
    exit 1
fi

echo "Starting FASTQ downloads..."

while read -r RUN; do

    if find "$OUT/fastq" -maxdepth 1 \
        -name "${RUN}*.fastq.gz" | grep -q .; then
        echo "Already downloaded: $RUN"
        continue
    fi

    echo "Downloading $RUN"

    SUCCESS=0

    for ATTEMPT in 1 2 3; do

        if prefetch "$RUN" \
            --output-directory "$OUT/sra" \
            --max-size 100G &&
           fasterq-dump "$OUT/sra/$RUN" \
            --split-files \
            --threads 4 \
            --outdir "$OUT/fastq" \
            --temp "$OUT/temp"; then

            SUCCESS=1
            break
        fi

        echo "Download retry $ATTEMPT: $RUN"
        sleep 15
    done

    if [[ "$SUCCESS" -eq 1 ]]; then

        find "$OUT/fastq" -maxdepth 1 \
            -name "${RUN}*.fastq" \
            -print0 | xargs -0 -r pigz -p 4

        rm -rf "$OUT/sra/$RUN"

        echo "Completed: $RUN"
    else
        echo "FAILED: $RUN; will retain it in the run list."
    fi

done < "$RUNFILE"

echo "RNA-seq download process finished."
echo "FASTQ directory: $OUT/fastq"

find "$OUT/fastq" -maxdepth 1 -name '*.fastq.gz' | wc -l
