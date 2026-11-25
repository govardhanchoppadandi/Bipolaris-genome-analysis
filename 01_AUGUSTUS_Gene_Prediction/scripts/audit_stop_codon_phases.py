import csv
from pathlib import Path
from collections import defaultdict, Counter

base = Path("test_predictions")
models = [
    "bipolaris_nd90pr",
    "bipolaris_nd90pr_v2",
    "bipolaris_multi_v3_all21",
    "bipolaris_multi_v4_advanced",
    "bipolaris_multi_v5_clean",
    "bipolaris_multi_v5_deduplicated",
]

stop_file = base / "internal_stop_genomic_coordinates.tsv"

with stop_file.open() as f:
    stops = list(csv.DictReader(f, delimiter="\t"))

stops_by_model = defaultdict(list)
for row in stops:
    stops_by_model[row["Model"]].append(row)

all_details = []

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
                "phase": int(fields[7]),
            })

    phase_mismatches = 0
    phase_checked = 0
    stop_patterns = Counter()

    for transcript, segments in features.items():
        strand = segments[0]["strand"]

        segments.sort(
            key=lambda x: x["start"],
            reverse=(strand == "-")
        )

        cumulative_length = 0

        for i, seg in enumerate(segments):
            expected = (
                (3 - cumulative_length % 3) % 3
                if i > 0 else seg["phase"]
            )

            if i > 0:
                phase_checked += 1
                if seg["phase"] != expected:
                    phase_mismatches += 1

            seg["expected_phase"] = expected
            seg["segment_index"] = i
            cumulative_length += seg["end"] - seg["start"] + 1

    for row in stops_by_model[model]:
        transcript = row["Sequence_ID"]
        segments = features.get(transcript, [])

        if not segments:
            continue

        strand = segments[0]["strand"]
        segments.sort(
            key=lambda x: x["start"],
            reverse=(strand == "-")
        )

        coords = [int(x) for x in row["Genomic_Coordinates"].split(",")]
        segment_hits = []

        for coord in coords:
            hit = None
            for i, seg in enumerate(segments):
                if seg["start"] <= coord <= seg["end"]:
                    hit = i
                    break
            segment_hits.append(hit)

        pattern = tuple(segment_hits)
        stop_patterns[pattern] += 1

        phase_info = []
        for idx in sorted(set(x for x in segment_hits if x is not None)):
            seg = segments[idx]
            phase_info.append(
                f"exon{idx+1}:phase={seg['phase']}:"
                f"expected={seg['expected_phase']}:"
                f"{seg['start']}-{seg['end']}"
            )

        all_details.append({
            "Model": model,
            "Sequence_ID": transcript,
            "Stop_Codon": row["Stop_Codon"],
            "Genomic_Coordinates": row["Genomic_Coordinates"],
            "Exon_Segment_Indices": ",".join(
                "NA" if x is None else str(x+1)
                for x in segment_hits
            ),
            "Exon_Distribution": ",".join(
                str(Counter(segment_hits).get(i, 0))
                for i in sorted(set(x for x in segment_hits if x is not None))
            ),
            "Segment_Phases": ";".join(phase_info),
        })

    print("\n" + "=" * 65)
    print(model)
    print("Transcripts with CDS segments:", len(features))
    print("Internal stop codons:", len(stops_by_model[model]))
    print("Non-first CDS segments phase-checked:", phase_checked)
    print("Phase mismatches:", phase_mismatches)
    print("Stop codon exon distributions:")

    for pattern, count in stop_patterns.most_common():
        print(" ", pattern, ":", count)

out = base / "stop_codon_phase_audit.tsv"

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=all_details[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(all_details)

print("\nDetailed report saved to:", out)
print("\nFirst 20 stop-codon phase records:")

for row in all_details[:20]:
    print(
        row["Model"],
        row["Sequence_ID"],
        row["Stop_Codon"],
        row["Exon_Segment_Indices"],
        row["Segment_Phases"]
    )
