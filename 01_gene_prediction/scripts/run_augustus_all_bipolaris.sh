#!/usr/bin/env bash
# Run AUGUSTUS predictions on accession directories directly inside each root.
# Usage:
#   bash run_augustus_all_bipolaris.sh "/path/ARI genomes" "/path/Bipolaris"
set -uo pipefail
export LC_ALL=C

AUGUSTUS_SPECIES="${AUGUSTUS_SPECIES:-bipolaris_nd90pr}"
ROOT1="${1:-/mnt/d/Bipolaris_genomes/ARI genomes}"
ROOT2="${2:-/mnt/d/Bipolaris_genomes/Bipolaris}"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: Activate the Conda environment first." >&2
  exit 1
fi
export AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH:-$CONDA_PREFIX/config}"

for program in augustus gffread find; do
  command -v "$program" >/dev/null 2>&1 || {
    echo "ERROR: Required command not found: $program" >&2
    exit 1
  }
done

if [[ ! -d "$AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES" ]]; then
  echo "ERROR: AUGUSTUS model not found: $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES" >&2
  exit 1
fi

for BASE_DIR in "$ROOT1" "$ROOT2"; do
  if [[ ! -d "$BASE_DIR" ]]; then
    echo "WARNING: Directory not found: $BASE_DIR" >&2
    continue
  fi

  REPORT="$BASE_DIR/gene_prediction_report.tsv"
  LOG="$BASE_DIR/augustus_prediction.log"
  printf 'Genome_Folder\tGenome_FASTA\tStatus\tReason\n' > "$REPORT"
  : > "$LOG"

  completed=0
  skipped=0
  failed=0
  missing=0
  total=0

  echo
  echo "========================================"
  echo "Directory: $BASE_DIR"
  echo "Model:     $AUGUSTUS_SPECIES"
  echo "========================================"

  while IFS= read -r -d '' asm; do
    ((total+=1))
    genome_name="$(basename "$asm")"

    # Skip only if all three outputs are present and non-empty.
    if [[ -s "$asm/genes.gff3" && -s "$asm/proteins.faa" && -s "$asm/cds.fna" ]]; then
      echo "SKIPPED: $genome_name (complete output set exists)"
      printf '%s\t\tSKIPPED\tComplete output set already exists\n' "$asm" >> "$REPORT"
      ((skipped+=1))
      continue
    fi

    fasta=""
    shopt -s nullglob
    candidates=( "$asm"/*.fna "$asm"/*.FNA "$asm"/*.fasta "$asm"/*.FASTA "$asm"/*.fa "$asm"/*.FA )
    shopt -u nullglob
    for f in "${candidates[@]}"; do
      if [[ -f "$f" ]]; then
        fasta="$f"
        break
      fi
    done

    if [[ -z "$fasta" ]]; then
      echo "NO FASTA: $genome_name"
      printf '%s\t\tNO_FASTA\tNo FASTA file directly inside accession folder\n' "$asm" >> "$REPORT"
      ((missing+=1))
      continue
    fi

    gff="$asm/genes.gff3"
    prot="$asm/proteins.faa"
    cds="$asm/cds.fna"
    echo "Running: $genome_name"
    echo "Genome:  $fasta"

    # Write to temporary files first so a failed run does not leave apparently complete outputs.
    tmpgff="$gff.tmp.$$"
    tmpprot="$prot.tmp.$$"
    tmpcds="$cds.tmp.$$"
    if augustus --species="$AUGUSTUS_SPECIES" --gff3=on --UTR=off "$fasta" \
        > "$tmpgff" 2>> "$LOG"; then
      if gffread "$tmpgff" -g "$fasta" -y "$tmpprot" -x "$tmpcds" 2>> "$LOG" &&
         [[ -s "$tmpgff" && -s "$tmpprot" && -s "$tmpcds" ]]; then
        mv "$tmpgff" "$gff"
        mv "$tmpprot" "$prot"
        mv "$tmpcds" "$cds"
        echo "COMPLETED: $genome_name"
        printf '%s\t%s\tCOMPLETED\tAUGUSTUS prediction and gffread extraction succeeded\n' \
          "$asm" "$fasta" >> "$REPORT"
        ((completed+=1))
      else
        rm -f "$tmpgff" "$tmpprot" "$tmpcds"
        echo "FAILED: gffread or empty output for $genome_name"
        printf '%s\t%s\tFAILED\tgffread failed or output empty; see log\n' \
          "$asm" "$fasta" >> "$REPORT"
        ((failed+=1))
      fi
    else
      rm -f "$tmpgff" "$tmpprot" "$tmpcds"
      echo "FAILED: AUGUSTUS for $genome_name"
      printf '%s\t%s\tFAILED\tAUGUSTUS failed; see log\n' "$asm" "$fasta" >> "$REPORT"
      ((failed+=1))
    fi
  done < <(find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

  echo
  echo "============== SUMMARY =============="
  echo "Directory:              $BASE_DIR"
  echo "Genome folders checked: $total"
  echo "Completed:              $completed"
  echo "Skipped:                $skipped"
  echo "Failed:                 $failed"
  echo "No genome FASTA:        $missing"
  echo "Report: $REPORT"
  echo "Log:    $LOG"
  echo "======================================"
done
