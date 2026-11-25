from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import csv

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub" / "test_predictions"

models = [
    "bipolaris_nd90pr",
    "bipolaris_nd90pr_v2",
    "bipolaris_multi_v3_all21",
    "bipolaris_multi_v4_advanced",
    "bipolaris_multi_v5_clean",
    "bipolaris_multi_v5_deduplicated",
]

output = base / "translation_mismatch_details.tsv"
rows = []

for model in models:
    proteins = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(base / "proteins" / f"{model}.faa", "fasta")
    }

    cds_records = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(base / "cds" / f"{model}.fna", "fasta")
    }

    for seq_id in sorted(set(proteins) & set(cds_records)):
        cds = cds_records[seq_id]
        protein = proteins[seq_id]

        if len(cds) % 3 != 0:
            rows.append({
                "Model": model,
                "Sequence_ID": seq_id,
                "Issue": "CDS_length_not_multiple_of_3",
                "CDS_length": len(cds),
                "Protein_length": len(protein),
                "Stop_codon_positions": "",
                "First_difference": ""
            })
            continue

        translated = str(Seq(cds).translate()).rstrip("*")
        protein_clean = protein.rstrip("*")

        stop_positions = [
            str(i + 1)
            for i, aa in enumerate(str(Seq(cds).translate()))
            if aa == "*"
        ]

        if translated != protein_clean:
            first_diff = next(
                (
                    str(i + 1)
                    for i, (a, b) in enumerate(zip(translated, protein_clean))
                    if a != b
                ),
                "length_difference"
            )

            rows.append({
                "Model": model,
                "Sequence_ID": seq_id,
                "Issue": "Translation_mismatch",
                "CDS_length": len(cds),
                "Protein_length": len(protein),
                "Stop_codon_positions": ",".join(stop_positions),
                "First_difference": first_diff
            })

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "Model", "Sequence_ID", "Issue", "CDS_length",
            "Protein_length", "Stop_codon_positions",
            "First_difference"
        ],
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("Diagnostic report:", output)
print("Total affected records:", len(rows))

for row in rows[:20]:
    print(row)
