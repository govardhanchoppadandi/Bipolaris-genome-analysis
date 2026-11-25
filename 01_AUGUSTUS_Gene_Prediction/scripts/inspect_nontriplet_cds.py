from Bio import SeqIO
from pathlib import Path

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub" / "test_predictions"

targets = {
    "bipolaris_multi_v3_all21": ["g5734.t1", "g9926.t1"],
    "bipolaris_multi_v4_advanced": ["g5586.t1"],
    "bipolaris_multi_v5_clean": ["g5713.t1", "g9907.t1"],
}

for model, ids in targets.items():
    print("\n" + "=" * 65)
    print("MODEL:", model)

    cds_file = base / "cds" / f"{model}.fna"

    for record in SeqIO.parse(cds_file, "fasta"):
        if record.id not in ids:
            continue

        seq = str(record.seq).upper()

        print("\nSequence:", record.id)
        print("CDS length:", len(seq))
        print("Length remainder:", len(seq) % 3)
        print("First 30 nt:", seq[:30])
        print("Last 60 nt:", seq[-60:])
        print("Trailing incomplete bases:", seq[-(len(seq) % 3):])

        complete_length = len(seq) - (len(seq) % 3)
        print("Complete codons:", complete_length // 3)

        print("Trailing complete codon:",
              seq[complete_length-3:complete_length]
              if complete_length >= 3 else "NA")
