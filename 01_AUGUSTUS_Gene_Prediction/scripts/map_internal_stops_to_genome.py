from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from collections import defaultdict
import csv

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub"
pred = base / "test_predictions"

genome_file = Path(
    "/mnt/d/Bipolaris_genomes/Bipolaris/GCA_000527765.2/"
    "GCA_000527765.2_ANU_Cvic_canu1.6_pilon5_genomic.fna"
)

models = [
    "bipolaris_nd90pr",
    "bipolaris_nd90pr_v2",
    "bipolaris_multi_v3_all21",
    "bipolaris_multi_v4_advanced",
    "bipolaris_multi_v5_clean",
    "bipolaris_multi_v5_deduplicated",
]

genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
results = []

for model in models:
    print("\n" + "=" * 65)
    print("MODEL:", model)

    cds_file = pred / "cds" / f"{model}.fna"
    gff_file = pred / f"{model}.gff3"

    cds_records = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(cds_file, "fasta")
    }

    features = defaultdict(list)

    with gff_file.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            f = line.rstrip().split("\t")
            if len(f) != 9 or f[2] != "CDS":
                continue

            attrs = dict(
                item.split("=", 1)
                for item in f[8].split(";")
                if "=" in item
            )

            parent = attrs.get("Parent")
            if parent:
                features[parent].append(f)

    ambiguous_records = set()
    mismatch_records = set()
    model_stops = 0

    for seq_id, cds in cds_records.items():
        ambiguous_count = sum(b not in "ACGT" for b in cds)

        if ambiguous_count:
            ambiguous_records.add(seq_id)

        translated = str(Seq(cds).translate())

        stop_positions = [
            i for i, aa in enumerate(translated[:-1])
            if aa == "*"
        ]

        if not stop_positions:
            continue

        model_stops += len(stop_positions)

        cds_features = features.get(seq_id, [])

        if not cds_features:
            print("WARNING: No GFF3 CDS features for", seq_id)
            continue

        strand = cds_features[0][6]
        contig = cds_features[0][0]

        cds_features.sort(
            key=lambda x: int(x[3]),
            reverse=(strand == "-")
        )

        reconstructed = ""
        coordinate_map = []

        for f in cds_features:
            start = int(f[3])
            end = int(f[4])

            segment = genome[contig].seq[start-1:end]

            if strand == "-":
                reconstructed += str(segment.reverse_complement()).upper()
                coordinate_map.extend(range(end, start-1, -1))
            else:
                reconstructed += str(segment).upper()
                coordinate_map.extend(range(start, end+1))

        if reconstructed != cds:
            mismatch_records.add(seq_id)

        for pos in stop_positions:
            codon = cds[pos*3:pos*3+3]
            coords = coordinate_map[pos*3:pos*3+3]

            left = max(0, pos*3-15)
            right = min(len(cds), pos*3+18)
            context = cds[left:right]

            results.append({
                "Model": model,
                "Sequence_ID": seq_id,
                "Contig": contig,
                "Strand": strand,
                "CDS_Codon_Position": pos+1,
                "Stop_Codon": codon,
                "Genomic_Coordinates": ",".join(map(str, coords)),
                "CDS_Ambiguous_Bases": ambiguous_count,
                "GFF3_Reconstruction_Matches": reconstructed == cds,
                "CDS_Context": context,
            })

    print("Affected CDS records:", len({
        r["Sequence_ID"] for r in results if r["Model"] == model
    }))
    print("Internal stop codons:", model_stops)
    print("CDS records with ambiguous bases:", len(ambiguous_records))
    print("CDS records not matching GFF3 reconstruction:",
          len(mismatch_records))

out = pred / "internal_stop_genomic_coordinates.tsv"

with out.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=results[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("\nOutput:", out)
print("\nFirst 20 mapped internal stops:")

for row in results[:20]:
    print(
        row["Model"],
        row["Sequence_ID"],
        row["Stop_Codon"],
        row["Genomic_Coordinates"],
        "GFF3 match:", row["GFF3_Reconstruction_Matches"]
    )
