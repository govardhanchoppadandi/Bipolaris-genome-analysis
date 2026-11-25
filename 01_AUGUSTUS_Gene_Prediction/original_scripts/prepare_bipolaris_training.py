from collections import defaultdict

REF = "/mnt/d/Bipolaris_genomes/Bipolaris/GCF_000338995.1"
GFF = REF + "/genomic.gff"
OUT = REF + "/genomic.complete_training.gff3"

genes = {}
transcripts = {}
transcript_parent = {}
cds_by_transcript = defaultdict(list)
records = []

def parse_attributes(s):
    result = {}
    for item in s.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            result[k] = v
    return result

with open(GFF) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue

        feature = fields[2]
        attr = parse_attributes(fields[8])
        records.append((fields, attr))

        if feature == "gene":
            genes[attr.get("ID")] = attr

        elif feature == "mRNA":
            tid = attr.get("ID")
            transcripts[tid] = attr
            transcript_parent[tid] = attr.get("Parent")

        elif feature == "CDS":
            for tid in attr.get("Parent", "").split(","):
                cds_by_transcript[tid].append(attr)

keep_genes = set()
keep_transcripts = set()

for tid, transcript in transcripts.items():
    gid = transcript_parent.get(tid)
    gene = genes.get(gid)
    cds = cds_by_transcript.get(tid, [])

    if not gene or not cds:
        continue

    partial = (
        gene.get("partial") == "true"
        or transcript.get("partial") == "true"
        or any(x.get("partial") == "true" for x in cds)
    )

    if not partial:
        keep_genes.add(gid)
        keep_transcripts.add(tid)

with open(OUT, "w") as out:
    out.write("##gff-version 3\n")

    for fields, attr in records:
        feature = fields[2]
        feature_id = attr.get("ID")
        parents = set(attr.get("Parent", "").split(","))

        if feature == "gene" and feature_id in keep_genes:
            out.write("\t".join(fields) + "\n")

        elif feature == "mRNA" and feature_id in keep_transcripts:
            out.write("\t".join(fields) + "\n")

        elif parents.intersection(keep_genes | keep_transcripts):
            out.write("\t".join(fields) + "\n")

print("Filtered GFF:", OUT)
print("Selected genes:", len(keep_genes))
print("Selected transcripts:", len(keep_transcripts))
