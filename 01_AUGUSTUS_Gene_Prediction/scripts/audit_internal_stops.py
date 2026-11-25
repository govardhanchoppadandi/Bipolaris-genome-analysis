from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import csv

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub" / "test_predictions"
out = base / "internal_stop_details.tsv"

models = [
    "bipolaris_nd90pr",
    "bipolaris_nd90pr_v2",
    "bipolaris_multi_v3_all21",
    "bipolaris_multi_v4_advanced",
    "bipolaris_multi_v5_clean",
    "bipolaris_multi_v5_deduplicated",
]

rows = []

for model in models:
    cds_file = base / "cds" / f"{model}.fna"

    for record in SeqIO.parse(cds_file, "fasta"):
        cds = str(record.seq).upper()

        if len(cds) % 3 != 0:
            continue

        translated = str(Seq(cds).translate())
        internal_positions = [
            i for i, aa in enumerate(translated[:-1])
            if aa == "*"
        ]

        for pos in internal_positions:
            codon = cds[pos * 3:pos * 3 + 3]

            rows.append({
                "Model": model,
                "Sequence_ID": record.id,
                "CDS_Length": len(cds),
                "Stop_Codon_Position": pos + 1,
                "Stop_Codon": codon,
                "Ambiguous_Bases_In_CDS": sum(
                    base not in "ACGT" for base in cds
                ),
                "Ambiguous_Bases_In_Stop_Codon": sum(
                    base not in "ACGT" for base in codon
                ),
                "Protein_Context": translated[
                    max(0, pos - 5):pos + 6
                ],
            })

with out.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=rows[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("Output:", out)
print("Total internal stop codons:", len(rows))

for model in models:
    selected = [r for r in rows if r["Model"] == model]
    affected = len(set(r["Sequence_ID"] for r in selected))
    print(model, "affected sequences:", affected,
          "internal stops:", len(selected))
