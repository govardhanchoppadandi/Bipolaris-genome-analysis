from Bio import SeqIO
from pathlib import Path
from collections import Counter, defaultdict
import csv
import statistics

BASE = Path.home() / "Bipolaris_AUGUSTUS_GitHub"
PRED = BASE / "test_predictions"

MODELS = [
    "bipolaris_nd90pr",
    "bipolaris_nd90pr_v2",
    "bipolaris_multi_v3_all21",
    "bipolaris_multi_v4_advanced",
    "bipolaris_multi_v5_clean",
    "bipolaris_multi_v5_deduplicated",
]

protein_summary = []
gene_summary = []
gene_details = []

for model in MODELS:
    protein_file = PRED / "proteins" / f"{model}.faa"
    cds_file = PRED / "cds" / f"{model}.fna"
    gff_file = PRED / f"{model}.gff3"

    proteins = list(SeqIO.parse(protein_file, "fasta"))
    cds_records = list(SeqIO.parse(cds_file, "fasta"))

    lengths = [len(r.seq) for r in proteins]

    # Parse GFF3 gene and CDS features
    genes = {}
    transcripts = {}
    cds_by_parent = defaultdict(list)
    start_codons = set()
    stop_codons = set()

    with gff_file.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            f = line.rstrip().split("\t")
            if len(f) != 9:
                continue

            feature = f[2]
            attrs = dict(
                item.split("=", 1)
                for item in f[8].split(";")
                if "=" in item
            )

            feature_id = attrs.get("ID")
            parent = attrs.get("Parent")

            if feature == "gene" and feature_id:
                genes[feature_id] = {
                    "contig": f[0],
                    "start": int(f[3]),
                    "end": int(f[4]),
                    "strand": f[6],
                }

            elif feature in ("mRNA", "transcript") and feature_id:
                transcripts[feature_id] = {
                    "contig": f[0],
                    "start": int(f[3]),
                    "end": int(f[4]),
                    "strand": f[6],
                    "parent": parent,
                }

            elif feature == "CDS" and parent:
                for p in parent.split(","):
                    cds_by_parent[p].append({
                        "start": int(f[3]),
                        "end": int(f[4]),
                        "phase": f[7],
                    })

            elif feature == "start_codon" and parent:
                start_codons.update(parent.split(","))

            elif feature == "stop_codon" and parent:
                stop_codons.update(parent.split(","))

    exon_counts = []
    cds_lengths = []
    missing_start = 0
    missing_stop = 0

    for transcript, info in transcripts.items():
        segments = cds_by_parent.get(transcript, [])
        if not segments:
            continue

        exon_count = len(segments)
        cds_length = sum(
            x["end"] - x["start"] + 1 for x in segments
        )

        exon_counts.append(exon_count)
        cds_lengths.append(cds_length)

        if transcript not in start_codons:
            missing_start += 1

        if transcript not in stop_codons:
            missing_stop += 1

        gene_details.append({
            "Model": model,
            "Transcript_ID": transcript,
            "Contig": info["contig"],
            "Strand": info["strand"],
            "Transcript_Length": info["end"] - info["start"] + 1,
            "CDS_Length": cds_length,
            "Exon_Count": exon_count,
            "Start_Codon_Feature": transcript in start_codons,
            "Stop_Codon_Feature": transcript in stop_codons,
        })

    def percentile(values, p):
        if not values:
            return 0
        return statistics.quantiles(
            values, n=100, method="inclusive"
        )[p - 1]

    protein_summary.append({
        "Model": model,
        "Protein_Count": len(lengths),
        "Minimum_aa": min(lengths) if lengths else 0,
        "Q1_aa": percentile(lengths, 25),
        "Median_aa": statistics.median(lengths) if lengths else 0,
        "Mean_aa": round(statistics.mean(lengths), 2) if lengths else 0,
        "Q3_aa": percentile(lengths, 75),
        "Maximum_aa": max(lengths) if lengths else 0,
        "Proteins_lt50aa": sum(x < 50 for x in lengths),
        "Proteins_lt100aa": sum(x < 100 for x in lengths),
        "Proteins_lt150aa": sum(x < 150 for x in lengths),
        "Proteins_gt1000aa": sum(x > 1000 for x in lengths),
        "Proteins_gt2000aa": sum(x > 2000 for x in lengths),
        "Proteins_with_internal_stop_symbol": sum(
            "*" in str(r.seq) or "." in str(r.seq) for r in proteins
        ),
        "CDS_Count": len(cds_records),
    })

    gene_summary.append({
        "Model": model,
        "Transcripts_with_CDS": len(exon_counts),
        "Mean_Exons_per_Transcript": round(
            statistics.mean(exon_counts), 2
        ) if exon_counts else 0,
        "Median_Exons_per_Transcript": statistics.median(
            exon_counts
        ) if exon_counts else 0,
        "Single_Exon_Transcripts": sum(x == 1 for x in exon_counts),
        "Multi_Exon_Transcripts": sum(x > 1 for x in exon_counts),
        "Median_CDS_Length_bp": statistics.median(
            cds_lengths
        ) if cds_lengths else 0,
        "Mean_CDS_Length_bp": round(
            statistics.mean(cds_lengths), 2
        ) if cds_lengths else 0,
        "Missing_Start_Codon_Feature": missing_start,
        "Missing_Stop_Codon_Feature": missing_stop,
    })

def write_tsv(filename, rows):
    path = PRED / filename
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=rows[0].keys(),
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    print("Saved:", path)

write_tsv("protein_length_summary.tsv", protein_summary)
write_tsv("gene_structure_summary.tsv", gene_summary)
write_tsv("gene_structure_details.tsv", gene_details)

print("\nPROTEIN LENGTH SUMMARY")
for row in protein_summary:
    print(row)

print("\nGENE STRUCTURE SUMMARY")
for row in gene_summary:
    print(row)
