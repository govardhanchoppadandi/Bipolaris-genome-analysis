#!/usr/bin/env python3
"""Copy one genomic FASTA and proteins.faa from each accession directory."""

from pathlib import Path
import csv
import shutil
import argparse

# Edit these defaults. For WSL, use Path("/mnt/d/Bipolaris_genomes/Bipolaris").
SOURCE_DIR = Path(r"D:\Bipolaris_genomes\Bipolaris")
BUSCO_DIR = SOURCE_DIR / "BUSCO"
GENOME_EXTENSIONS = {".fasta", ".fa", ".fna", ".fas"}


def inside(path, parent):
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=SOURCE_DIR)
    ap.add_argument("--busco-dir", type=Path, default=None)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    source = args.source.resolve()
    project = (args.busco_dir or source / "BUSCO").resolve()
    genome_out = project / "Genome input"
    protein_out = project / "Protein input"
    genome_out.mkdir(parents=True, exist_ok=True)
    protein_out.mkdir(parents=True, exist_ok=True)

    rows, genomes, proteins = [], set(), set()
    accession_dirs = sorted(
        p for p in source.iterdir()
        if p.is_dir() and not inside(p, project)
    )

    for folder in accession_dirs:
        accession = folder.name
        found = [
            p for p in folder.rglob("*")
            if p.is_file() and not inside(p, project)
        ]
        fasta = [p for p in found if p.suffix.lower() in GENOME_EXTENSIONS]
        preferred = [
            p for p in fasta
            if "genomic" in p.name.lower()
            and not p.name.lower().startswith(("cds", "rna"))
        ]
        candidates = preferred or [
            p for p in fasta
            if not p.name.lower().startswith(("cds", "rna"))
            and "protein" not in p.name.lower()
        ]
        protein_candidates = [p for p in found if p.name.lower() == "proteins.faa"]

        gs, ps, gsrc, psrc = "MISSING", "MISSING", "", ""
        if len(candidates) == 1:
            src = candidates[0]
            dst = genome_out / f"{accession}.fasta"
            gsrc = str(src)
            if dst.exists() and not args.overwrite:
                gs = "EXISTS_NOT_OVERWRITTEN"
            else:
                shutil.copy2(src, dst)
                gs = "COPIED"
            genomes.add(accession)
        elif len(candidates) > 1:
            gs = "MULTIPLE_CANDIDATES: " + " | ".join(map(str, candidates))

        if len(protein_candidates) == 1:
            src = protein_candidates[0]
            dst = protein_out / f"{accession}.faa"
            psrc = str(src)
            if dst.exists() and not args.overwrite:
                ps = "EXISTS_NOT_OVERWRITTEN"
            else:
                shutil.copy2(src, dst)
                ps = "COPIED"
            proteins.add(accession)
        elif len(protein_candidates) > 1:
            ps = "MULTIPLE proteins.faa FILES"

        rows.append({
            "Accession": accession, "Genome_status": gs, "Genome_source": gsrc,
            "Protein_status": ps, "Protein_source": psrc
        })

    project.mkdir(parents=True, exist_ok=True)
    report = project / "input_preparation_report.tsv"
    with report.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else [
            "Accession", "Genome_status", "Genome_source",
            "Protein_status", "Protein_source"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print("Accession folders:", len(rows))
    print("Genome accessions:", len(genomes))
    print("Protein accessions:", len(proteins))
    print("Matching accessions:", genomes == proteins)
    print("Only genome:", sorted(genomes - proteins))
    print("Only protein:", sorted(proteins - genomes))
    print("Review MULTIPLE_CANDIDATES in:", report)


if __name__ == "__main__":
    main()
