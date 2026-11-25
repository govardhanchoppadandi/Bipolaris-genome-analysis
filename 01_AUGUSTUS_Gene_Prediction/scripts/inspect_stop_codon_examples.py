from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub" / "test_predictions"

models = [
    "bipolaris_nd90pr",
    "bipolaris_nd90pr_v2",
    "bipolaris_multi_v3_all21",
    "bipolaris_multi_v4_advanced",
    "bipolaris_multi_v5_clean",
    "bipolaris_multi_v5_deduplicated",
]

for model in models:
    print("\n" + "=" * 70)
    print("MODEL:", model)

    proteins = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(base / "proteins" / f"{model}.faa", "fasta")
    }

    cds_records = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(base / "cds" / f"{model}.fna", "fasta")
    }

    shown = 0

    for seq_id, cds in cds_records.items():
        if len(cds) % 3 != 0:
            print("NON_MULTIPLE_OF_3:", seq_id, "CDS length:", len(cds))
            continue

        translated = str(Seq(cds).translate())
        protein = proteins.get(seq_id, "")

        internal_stops = [
            i for i, aa in enumerate(translated[:-1]) if aa == "*"
        ]

        if internal_stops:
            print("\nSequence:", seq_id)
            print("CDS length:", len(cds))
            print("Protein length:", len(protein))

            for pos in internal_stops[:3]:
                print("Stop position:", pos + 1)
                print("CDS codon:", cds[pos * 3:pos * 3 + 3])
                print("Protein amino acid:", protein[pos] if pos < len(protein) else "MISSING")
                print("Translated amino acid:", translated[pos])
                print("Protein context:", protein[max(0, pos-5):pos+6])

            shown += 1
            if shown >= 3:
                break
