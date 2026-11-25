from pathlib import Path

base = Path.home() / "Bipolaris_AUGUSTUS_GitHub" / "test_predictions"

targets = {
    "bipolaris_multi_v3_all21": ["g5734.t1", "g9926.t1"],
    "bipolaris_multi_v4_advanced": ["g5586.t1"],
    "bipolaris_multi_v5_clean": ["g5713.t1", "g9907.t1"],
}

for model, ids in targets.items():
    gff = base / f"{model}.gff3"

    print("\n" + "=" * 70)
    print("MODEL:", model)

    with gff.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            feature = fields[2]
            attributes = fields[8]

            for seq_id in ids:
                if seq_id in attributes:
                    print(line.rstrip())
