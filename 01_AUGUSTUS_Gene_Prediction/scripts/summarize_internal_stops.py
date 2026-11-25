import csv
from collections import Counter, defaultdict
from pathlib import Path

file = Path("test_predictions/internal_stop_details.tsv")

with file.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

print("TOTAL INTERNAL STOPS:", len(rows))

print("\nSTOP CODON DISTRIBUTION")
codons = Counter(row["Stop_Codon"] for row in rows)
for codon, count in sorted(codons.items()):
    print(codon, count)

print("\nAMBIGUOUS BASES")
ambiguous = defaultdict(set)

for row in rows:
    if int(row["Ambiguous_Bases_In_CDS"]) > 0:
        ambiguous[row["Model"]].add(row["Sequence_ID"])

for model in sorted(set(row["Model"] for row in rows)):
    affected = {
        row["Sequence_ID"] for row in rows
        if row["Model"] == model
    }
    print(
        model,
        "affected CDS:", len(affected),
        "with ambiguous bases:", len(ambiguous[model]),
        "without ambiguous bases:", len(affected - ambiguous[model])
    )

print("\nSTOP CODONS BY MODEL")
for model in sorted(set(row["Model"] for row in rows)):
    selected = [r for r in rows if r["Model"] == model]
    counts = Counter(r["Stop_Codon"] for r in selected)
    print(model, dict(counts))

print("\nFIRST 15 INTERNAL STOP RECORDS")
for row in rows[:15]:
    print(
        row["Model"],
        row["Sequence_ID"],
        "position:", row["Stop_Codon_Position"],
        "codon:", row["Stop_Codon"],
        "ambiguous bases:", row["Ambiguous_Bases_In_CDS"]
    )
