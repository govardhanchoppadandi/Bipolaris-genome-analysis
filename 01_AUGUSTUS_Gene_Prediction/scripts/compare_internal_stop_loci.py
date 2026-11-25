import csv
from collections import defaultdict, Counter
from pathlib import Path

source = Path("test_predictions/internal_stop_genomic_coordinates.tsv")

with source.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

# Use contig, strand, stop codon, and sorted genomic coordinates
# to identify matching genomic stop-codon loci across models.
loci = defaultdict(list)

for row in rows:
    coords = tuple(sorted(
        int(x) for x in row["Genomic_Coordinates"].split(",")
    ))

    key = (
        row["Contig"],
        row["Strand"],
        row["Stop_Codon"],
        coords
    )

    loci[key].append(row)

summary = Counter()
output = []

for key, records in loci.items():
    contig, strand, codon, coords = key
    models = sorted(set(r["Model"] for r in records))

    summary[len(models)] += 1

    output.append({
        "Contig": contig,
        "Strand": strand,
        "Stop_Codon": codon,
        "Genomic_Coordinates": ",".join(map(str, coords)),
        "Number_Of_Models": len(models),
        "Models": ";".join(models),
        "Transcript_IDs": ";".join(
            f'{r["Model"]}:{r["Sequence_ID"]}' for r in records
        ),
    })

output.sort(
    key=lambda r: (
        -r["Number_Of_Models"],
        r["Contig"],
        r["Genomic_Coordinates"]
    )
)

print("INTERNAL STOP LOCUS RECURRENCE")
print("-" * 60)

for n in sorted(summary, reverse=True):
    print(f"Found in {n} model(s): {summary[n]} genomic loci")

print("\nTotal stop records:", len(rows))
print("Unique genomic stop loci:", len(loci))

print("\nLOCUS SHARING ACROSS MODELS")
print("-" * 60)

for row in output[:30]:
    print(
        row["Number_Of_Models"],
        row["Stop_Codon"],
        row["Contig"],
        row["Genomic_Coordinates"],
        row["Models"]
    )

out = Path("test_predictions/internal_stop_locus_recurrence.tsv")

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=output[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(output)

print("\nSaved:", out)
