import csv
from pathlib import Path
from collections import Counter, defaultdict
from Bio import SeqIO

base = Path("test_predictions")

genome_file = Path(
    "/mnt/d/Bipolaris_genomes/Bipolaris/GCA_000527765.2/"
    "GCA_000527765.2_ANU_Cvic_canu1.6_pilon5_genomic.fna"
)

genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

with (base / "internal_stop_genomic_coordinates.tsv").open() as f:
    stops = list(csv.DictReader(f, delimiter="\t"))

models = sorted(set(row["Model"] for row in stops))
results = []

for model in models:
    gff = base / f"{model}.gff3"
    features = defaultdict(list)

    with gff.open() as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue

            attrs = dict(
                item.split("=", 1)
                for item in fields[8].split(";")
                if "=" in item
            )

            parent = attrs.get("Parent")
            if not parent:
                continue

            features[parent].append({
                "contig": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "phase": fields[7],
            })

    counts = Counter()

    for row in [r for r in stops if r["Model"] == model]:
        transcript = row["Sequence_ID"]
        segments = features.get(transcript, [])

        if not segments:
            continue

        strand = segments[0]["strand"]
        contig = segments[0]["contig"]

        segments.sort(
            key=lambda x: x["start"],
            reverse=(strand == "-")
        )

        coords = [
            int(x) for x in row["Genomic_Coordinates"].split(",")
        ]

        segment_hits = []

        for coord in coords:
            hit = next(
                (
                    i for i, seg in enumerate(segments)
                    if seg["start"] <= coord <= seg["end"]
                ),
                None
            )
            segment_hits.append(hit)

        # Find the CDS segment transition crossed by this stop codon.
        transitions = []
        for a, b in zip(segment_hits, segment_hits[1:]):
            if a is not None and b is not None and a != b:
                transitions.append((a, b))

        if not transitions:
            continue

        for a, b in transitions:
            upstream = segments[a]
            downstream = segments[b]

            if strand == "+":
                intron_start = upstream["end"] + 1
                intron_end = downstream["start"] - 1
            else:
                intron_start = downstream["end"] + 1
                intron_end = upstream["start"] - 1

            if intron_start > intron_end:
                motif = "INVALID_INTRON_INTERVAL"
                donor = acceptor = ""
            else:
                intron = genome[contig].seq[
                    intron_start - 1:intron_end
                ]

                if strand == "-":
                    intron = intron.reverse_complement()

                intron = str(intron).upper()
                donor = intron[:2]
                acceptor = intron[-2:]
                motif = donor + "-" + acceptor

            if motif in ("GT-AG", "GC-AG", "AT-AC"):
                classification = "recognized_canonical"
            else:
                classification = "other_or_noncanonical"

            counts[classification] += 1

            results.append({
                "Model": model,
                "Sequence_ID": transcript,
                "Stop_Codon": row["Stop_Codon"],
                "Stop_Genomic_Coordinates": row["Genomic_Coordinates"],
                "Strand": strand,
                "Intron_Start": intron_start,
                "Intron_End": intron_end,
                "Donor": donor,
                "Acceptor": acceptor,
                "Splice_Motif": motif,
                "Classification": classification,
            })

    print("\n", model)
    print("Stop-codon junctions examined:", sum(counts.values()))
    print("Recognized canonical:", counts["recognized_canonical"])
    print("Other/noncanonical:", counts["other_or_noncanonical"])

    motif_counts = Counter(
        r["Splice_Motif"]
        for r in results
        if r["Model"] == model
    )
    print("Motifs:", dict(motif_counts))

out = base / "stop_junction_splice_motifs.tsv"

if results:
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=results[0].keys(),
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(results)

print("\nDetailed report:", out)
print("\nFirst 25 junctions:")

for row in results[:25]:
    print(
        row["Model"],
        row["Sequence_ID"],
        row["Stop_Codon"],
        row["Splice_Motif"],
        row["Classification"]
    )
