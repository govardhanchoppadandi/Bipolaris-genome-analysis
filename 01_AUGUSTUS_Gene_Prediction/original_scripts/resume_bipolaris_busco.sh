# ============================================================
# BIPOLARIS BUSCO 6.1.0 - RESUME PIPELINE
# Lineage: pleosporales_odb12
# Genome: 74 FASTA files
# Protein: 74 FAA files
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate busco_latest

set -Eeuo pipefail

BASE="/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO"
LINEAGE="pleosporales_odb12"
DOWNLOAD="$HOME/busco_downloads"

GENOME_INPUT="$BASE/Genome input"
PROTEIN_INPUT="$BASE/Protein input"

GENOME_OUT="$BASE/BUSCO_RESULTS/GENOME"
PROTEIN_OUT="$BASE/BUSCO_RESULTS/PROTEIN"
FAILED="$BASE/BUSCO_RESULTS/FAILED"
LOGS="$BASE/BUSCO_RESULTS/LOGS"
RESULTS="$BASE/results"

THREADS=8

echo "=================================================="
echo "       BIPOLARIS BUSCO RESUME PIPELINE"
echo "=================================================="

# 1. Verify BUSCO

busco --version
echo "BUSCO location: $(which busco)"

if [[ ! -d "$DOWNLOAD/lineages/$LINEAGE" ]]; then
    echo "ERROR: Pleosporales lineage not found."
    exit 1
fi

echo "Lineage verified: $LINEAGE"

# 2. Create output folders

mkdir -p "$GENOME_OUT" "$PROTEIN_OUT" \
         "$FAILED" "$LOGS" "$RESULTS"

# 3. Check input files

shopt -s nullglob

GENOMES=("$GENOME_INPUT"/*.fasta)
PROTEINS=("$PROTEIN_INPUT"/*.faa)

echo "Genome files: ${#GENOMES[@]}"
echo "Protein files: ${#PROTEINS[@]}"

if [[ ${#GENOMES[@]} -eq 0 || ${#PROTEINS[@]} -eq 0 ]]; then
    echo "ERROR: Input files missing."
    exit 1
fi

# Function: find and retain only summary files

clean_successful_output() {
    local OUT="$1"
    local NAME="$2"
    local TYPE="$3"

    local TEMP="$OUT/summary_temp"
    mkdir -p "$TEMP"

    local TXT JSON

    TXT=$(find "$OUT" -type f -name 'short_summary*.txt' \
          ! -path "$TEMP/*" -print -quit)

    JSON=$(find "$OUT" -type f -name 'short_summary*.json' \
           ! -path "$TEMP/*" -print -quit)

    if [[ -z "$TXT" || -z "$JSON" ]]; then
        rm -rf "$TEMP"
        return 1
    fi

    cp "$TXT" "$TEMP/"
    cp "$JSON" "$TEMP/"

    find "$OUT" -mindepth 1 -maxdepth 1 \
        ! -name summary_temp -exec rm -rf {} +

    mv "$TEMP"/* "$OUT/"
    rmdir "$TEMP"

    echo "SUCCESS: $TYPE $NAME"
    return 0
}

# Function: check whether both summaries exist

is_completed() {
    local OUT="$1"

    [[ -n "$(find "$OUT" -maxdepth 1 -type f \
        -name 'short_summary*.txt' -size +0c -print -quit 2>/dev/null)" ]] &&
    [[ -n "$(find "$OUT" -maxdepth 1 -type f \
        -name 'short_summary*.json' -size +0c -print -quit 2>/dev/null)" ]]
}

# ============================================================
# 4. RESUME GENOME ANALYSIS
# ============================================================

echo
echo "=================================================="
echo "          RESUMING GENOME ANALYSIS"
echo "=================================================="

G_DONE=0
G_RUN=0
G_FAIL=0

for FILE in "${GENOMES[@]}"; do

    NAME=$(basename "$FILE" .fasta)
    OUT="$GENOME_OUT/$NAME"
    LOG="$LOGS/GENOME_${NAME}.log"

    if [[ -d "$OUT" ]] && is_completed "$OUT"; then
        echo "ALREADY COMPLETED - SKIP: $NAME"
        ((G_DONE+=1))
        continue
    fi

    # Preserve incomplete/interrupted previous output
    if [[ -d "$OUT" ]]; then
        DEST="$FAILED/GENOME_${NAME}_INTERRUPTED_$(date +%Y%m%d_%H%M%S)"
        mv "$OUT" "$DEST"
        echo "Previous incomplete output preserved: $DEST"
    fi

    echo
    echo "RUNNING GENOME: $NAME"

    if busco \
        -i "$FILE" \
        -m genome \
        -l "$LINEAGE" \
        -o "$NAME" \
        --out_path "$GENOME_OUT" \
        -c "$THREADS" \
        --download_path "$DOWNLOAD" \
        > "$LOG" 2>&1; then

        if clean_successful_output "$OUT" "$NAME" "GENOME"; then
            ((G_RUN+=1))
        else
            DEST="$FAILED/GENOME_${NAME}_SUMMARY_MISSING_$(date +%Y%m%d_%H%M%S)"
            [[ -d "$OUT" ]] && mv "$OUT" "$DEST"
            echo "FAILED: Missing TXT or JSON: $NAME"
            ((G_FAIL+=1))
        fi

    else
        DEST="$FAILED/GENOME_${NAME}_FAILED_$(date +%Y%m%d_%H%M%S)"
        [[ -d "$OUT" ]] && mv "$OUT" "$DEST"
        echo "FAILED: $NAME"
        ((G_FAIL+=1))
    fi

done

echo
echo "GENOME STAGE FINISHED"
echo "Previously completed: $G_DONE"
echo "Newly completed: $G_RUN"
echo "Failed: $G_FAIL"

# ============================================================
# 5. RESUME PROTEIN ANALYSIS
# ============================================================

echo
echo "=================================================="
echo "          RESUMING PROTEIN ANALYSIS"
echo "=================================================="

P_DONE=0
P_RUN=0
P_FAIL=0

for FILE in "${PROTEINS[@]}"; do

    NAME=$(basename "$FILE" .faa)
    OUT="$PROTEIN_OUT/$NAME"
    LOG="$LOGS/PROTEIN_${NAME}.log"

    if [[ -d "$OUT" ]] && is_completed "$OUT"; then
        echo "ALREADY COMPLETED - SKIP: $NAME"
        ((P_DONE+=1))
        continue
    fi

    # Preserve incomplete/interrupted previous output
    if [[ -d "$OUT" ]]; then
        DEST="$FAILED/PROTEIN_${NAME}_INTERRUPTED_$(date +%Y%m%d_%H%M%S)"
        mv "$OUT" "$DEST"
        echo "Previous incomplete output preserved: $DEST"
    fi

    echo
    echo "RUNNING PROTEIN: $NAME"

    if busco \
        -i "$FILE" \
        -m proteins \
        -l "$LINEAGE" \
        -o "$NAME" \
        --out_path "$PROTEIN_OUT" \
        -c "$THREADS" \
        --download_path "$DOWNLOAD" \
        > "$LOG" 2>&1; then

        if clean_successful_output "$OUT" "$NAME" "PROTEIN"; then
            ((P_RUN+=1))
        else
            DEST="$FAILED/PROTEIN_${NAME}_SUMMARY_MISSING_$(date +%Y%m%d_%H%M%S)"
            [[ -d "$OUT" ]] && mv "$OUT" "$DEST"
            echo "FAILED: Missing TXT or JSON: $NAME"
            ((P_FAIL+=1))
        fi

    else
        DEST="$FAILED/PROTEIN_${NAME}_FAILED_$(date +%Y%m%d_%H%M%S)"
        [[ -d "$OUT" ]] && mv "$OUT" "$DEST"
        echo "FAILED: $NAME"
        ((P_FAIL+=1))
    fi

done

echo
echo "PROTEIN STAGE FINISHED"
echo "Previously completed: $P_DONE"
echo "Newly completed: $P_RUN"
echo "Failed: $P_FAIL"

# ============================================================
# 6. GENERATE COMBINED EXCEL AND TSV REPORTS
# ============================================================

echo
echo "=================================================="
echo "          GENERATING FINAL REPORTS"
echo "=================================================="

python - <<'PY'
import re
import csv
from pathlib import Path

BASE = Path("/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO")
RESULTS = BASE / "results"
RESULTS.mkdir(exist_ok=True)

def read_busco(folder):
    records = {}

    for directory in folder.iterdir():
        if not directory.is_dir():
            continue

        files = list(directory.glob("short_summary*.txt"))
        if not files:
            continue

        text = files[0].read_text(errors="ignore")

        match = re.search(
            r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
            text
        )

        if match:
            c, s, d, f, m, n = match.groups()

            records[directory.name] = {
                "Complete_%": float(c),
                "Single_copy_%": float(s),
                "Duplicated_%": float(d),
                "Fragmented_%": float(f),
                "Missing_%": float(m),
                "BUSCO_count": int(n)
            }

    return records

genome = read_busco(BASE / "BUSCO_RESULTS" / "GENOME")
protein = read_busco(BASE / "BUSCO_RESULTS" / "PROTEIN")

columns = [
    "Sample",
    "Genome_Complete_%",
    "Genome_Single_copy_%",
    "Genome_Duplicated_%",
    "Genome_Fragmented_%",
    "Genome_Missing_%",
    "Genome_BUSCO_count",
    "Protein_Complete_%",
    "Protein_Single_copy_%",
    "Protein_Duplicated_%",
    "Protein_Fragmented_%",
    "Protein_Missing_%",
    "Protein_BUSCO_count",
    "Classification"
]

rows = []

for sample in sorted(set(genome) | set(protein)):
    g = genome.get(sample, {})
    p = protein.get(sample, {})

    if not g or not p:
        classification = "INCOMPLETE"
    elif (
        g["Complete_%"] >= 95 and p["Complete_%"] >= 95
        and g["Fragmented_%"] <= 2 and p["Fragmented_%"] <= 2
    ):
        classification = "KEEP"
    elif (
        g["Complete_%"] < 90 or p["Complete_%"] < 90
        or g["Fragmented_%"] > 5 or p["Fragmented_%"] > 5
    ):
        classification = "REMOVE_FRAGMENTED"
    else:
        classification = "REVIEW"

    rows.append({
        "Sample": sample,
        **{"Genome_" + k: v for k, v in g.items()},
        **{"Protein_" + k: v for k, v in p.items()},
        "Classification": classification
    })

# TSV report
tsv = RESULTS / "Bipolaris_BUSCO_Genome_Protein_Summary.tsv"

with tsv.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("TSV report created.")

# Excel report
try:
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    ws = wb.active
    ws.title = "BUSCO Summary"

    ws.append(columns)

    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="1F4E78")

    for row in rows:
        ws.append([row.get(c, "") for c in columns])

    ws.freeze_panes = "B2"
    ws.auto_filter.ref = ws.dimensions

    for col in ws.columns:
        letter = get_column_letter(col[0].column)
        width = min(max(len(str(c.value or "")) for c in col) + 2, 30)
        ws.column_dimensions[letter].width = width

    wb.save(RESULTS / "Bipolaris_BUSCO_Genome_Protein_Summary.xlsx")
    print("Excel report created.")

except ImportError:
    print("WARNING: openpyxl unavailable; TSV created only.")

print("Genome summaries:", len(genome))
print("Protein summaries:", len(protein))
print("Combined samples:", len(rows))
print("Reports saved:", RESULTS)
PY

# ============================================================
# 7. FINAL STATUS
# ============================================================

echo
echo "=================================================="
echo "             RESUME PIPELINE REPORT"
echo "=================================================="

echo "Genome TXT summaries:"
find "$GENOME_OUT" -type f -name 'short_summary*.txt' | wc -l

echo "Genome JSON summaries:"
find "$GENOME_OUT" -type f -name 'short_summary*.json' | wc -l

echo "Protein TXT summaries:"
find "$PROTEIN_OUT" -type f -name 'short_summary*.txt' | wc -l

echo "Protein JSON summaries:"
find "$PROTEIN_OUT" -type f -name 'short_summary*.json' | wc -l

echo "Failed/interrupted folders:"
find "$FAILED" -mindepth 1 -maxdepth 1 -type d | wc -l

echo
echo "FINAL REPORTS:"
ls -lh "$RESULTS"

echo
echo "=================================================="
echo "          RESUME PIPELINE FINISHED"
echo "=================================================="
