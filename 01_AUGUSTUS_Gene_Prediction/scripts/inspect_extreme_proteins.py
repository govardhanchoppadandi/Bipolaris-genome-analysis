from Bio import SeqIO
from pathlib import Path
from collections import defaultdict
import csv

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub"
pred = base / "test_predictions"

models = [
    "bipolaris_nd90pr",
    "bipolaris_nd90pr_v2",
    "bipolaris_multi_v3_all21",
    "bipolaris_multi_v4_advanced",
    "bipolaris_multi_v5_clean",
    "bipolaris_multi_v5_deduplicated",
]

results = []

for model in models:
    protein_file = pred / "proteins" / f"{model}.faa"
    gff_file = pred / f"{model}.gff3"

    lengths = {
        r.id: len(r.seq)
        for r in SeqIO.parse(protein_file, "fasta")
    }

    cds = defaultdict(list)
    transcript_info = {}

    with gff_file.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            f = line.rstrip().split("\t")
            if len(f) != 9:
                continue

            attrs = dict(
                item.split("=", 1)
                for item in f[8].split(";")
                if "=" in item
            )

            feature = f[2]
            feature_id = attrs.get("ID")
            parent = attrs.get("Parent")

            if feature in ("mRNA", "transcript") and feature_id:
                transcript_info[feature_id] = {
                    "contig": f[0],
                    "start": int(f[3]),
                    "end": int(f[4]),
                    "strand": f[6],
                }

            elif feature == "CDS" and parent:
                for p in parent.split(","):
                    cds[p].append((int(f[3]), int(f[4])))

    print("\n" + "=" * 70)
    print(model)

    for transcript, length in sorted(
        lengths.items(), key=lambda x: x[1], reverse=True
    )[:10]:
        info = transcript_info.get(transcript, {})
        segments = cds.get(transcript, [])

        row = {
            "Model": model,
            "Transcript_ID": transcript,
            "Protein_Length_aa": length,
            "Contig": info.get("contig", ""),
            "Transcript_Start": info.get("start", ""),
            "Transcript_End": info.get("end", ""),
            "Strand": info.get("strand", ""),
            "Exon_Count": len(segments),
            "CDS_Length_bp": sum(e-s+1 for s, e in segments),
        }
        results.append(row)

        print(
            transcript,
            "length:", length,
            "contig:", row["Contig"],
            "strand:", row["Strand"],
            "exons:", row["Exon_Count"],
            "CDS bp:", row["CDS_Length_bp"],
            "coordinates:", row["Transcript_Start"],
            "-", row["Transcript_End"]
        )

out = pred / "extreme_protein_details.tsv"

with out.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=results[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("\nSaved:", out)
