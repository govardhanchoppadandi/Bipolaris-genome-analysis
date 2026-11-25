from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import csv

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub" / "test_predictions"
out = base / "sequence_integrity_audit_corrected.tsv"

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
    proteins = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(base / "proteins" / f"{model}.faa", "fasta")
    }

    cds_records = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(base / "cds" / f"{model}.fna", "fasta")
    }

    common = set(proteins) & set(cds_records)

    non_multiple = 0
    internal_stops = 0
    mismatches = 0

    for seq_id in common:
        cds = cds_records[seq_id]
        protein = proteins[seq_id].replace(".", "*")

        if len(cds) % 3 != 0:
            non_multiple += 1
            continue

        translated = str(Seq(cds).translate())

        internal = translated[:-1].count("*")
        internal_stops += internal

        translated = translated.rstrip("*")
        protein = protein.rstrip("*")

        if translated != protein:
            mismatches += 1

    rows.append({
        "Model": model,
        "Protein_records": len(proteins),
        "CDS_records": len(cds_records),
        "Matching_IDs": len(common),
        "CDS_length_not_multiple_of_3": non_multiple,
        "Internal_stop_codons": internal_stops,
        "Translation_mismatches_after_normalization": mismatches,
    })

with out.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=rows[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("Corrected audit:", out)

for row in rows:
    print(row)
