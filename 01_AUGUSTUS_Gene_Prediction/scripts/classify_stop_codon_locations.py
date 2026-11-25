import csv
from collections import defaultdict
from pathlib import Path

source = Path("test_predictions/internal_stop_genomic_coordinates.tsv")

with source.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

summary = defaultdict(lambda: {
    "total": 0,
    "within_exon": 0,
    "crosses_exon_junction": 0
})

junction_rows = []

for row in rows:
    coords = [int(x) for x in row["Genomic_Coordinates"].split(",")]
    strand = row["Strand"]

    expected_step = -1 if strand == "-" else 1

    crosses_junction = any(
        coords[i + 1] - coords[i] != expected_step
        for i in range(len(coords) - 1)
    )

    model = row["Model"]
    summary[model]["total"] += 1

    if crosses_junction:
        summary[model]["crosses_exon_junction"] += 1
        junction_rows.append(row)
    else:
        summary[model]["within_exon"] += 1

print("STOP CODON LOCATION SUMMARY")
print("-" * 75)

for model, values in sorted(summary.items()):
    print(
        model,
        "total:", values["total"],
        "within exon:", values["within_exon"],
        "crosses exon junction:", values["crosses_exon_junction"]
    )

print("\nSTOP CODONS CROSSING EXON JUNCTIONS")
print("-" * 75)

for row in junction_rows[:30]:
    print(
        row["Model"],
        row["Sequence_ID"],
        row["Stop_Codon"],
        row["Genomic_Coordinates"]
    )

out = Path("test_predictions/stop_codons_crossing_exon_junctions.tsv")

with out.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=rows[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(junction_rows)

print("\nJunction-crossing records:", len(junction_rows))
print("Saved:", out)
