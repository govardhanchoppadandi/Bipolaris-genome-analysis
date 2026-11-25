from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

genome = Path(
    "/mnt/d/Bipolaris_genomes/Bipolaris/"
    "GCA_000527765.2/"
    "GCA_000527765.2_ANU_Cvic_canu1.6_pilon5_genomic.fna"
)

gff_base = Path.home() / "Bipolaris_AUGUSTUS_GitHub" / "test_predictions"

targets = {
    "bipolaris_multi_v3_all21": ["g5734.t1", "g9926.t1"],
    "bipolaris_multi_v4_advanced": ["g5586.t1"],
    "bipolaris_multi_v5_clean": ["g5713.t1", "g9907.t1"],
}

sequences = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

for model, ids in targets.items():
    print("\n" + "=" * 75)
    print("MODEL:", model)

    gff = gff_base / f"{model}.gff3"

    for seq_id in ids:
        print("\nTARGET:", seq_id)

        features = []

        with gff.open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                fields = line.rstrip().split("\t")
                if len(fields) != 9:
                    continue

                attrs = fields[8]

                if (
                    f"ID={seq_id};" in attrs
                    or f"Parent={seq_id};" in attrs
                    or attrs.endswith(f"Parent={seq_id}")
                ):
                    features.append(fields)

        if not features:
            print("No features found")
            continue

        contig = features[0][0]
        seq = sequences[contig].seq
        print("Contig:", contig)
        print("Contig length:", len(seq))

        for f in features:
            print(
                "Feature:", f[2],
                "Start:", f[3],
                "End:", f[4],
                "Strand:", f[6],
                "Phase:", f[7],
                "Attributes:", f[8]
            )

        cds = [f for f in features if f[2] == "CDS"]

        if cds:
            low = min(int(f[3]) for f in cds)
            high = max(int(f[4]) for f in cds)

            print("CDS genomic span:", low, "-", high)
            print("Sequence around lower CDS boundary:")

            start = max(0, low - 31)
            end = min(len(seq), low + 30)

            region = seq[start:end]
            print("Forward:", region)
            print("Reverse complement:", region.reverse_complement())

        print("First 60 bases of contig:")
        print(seq[:60])
