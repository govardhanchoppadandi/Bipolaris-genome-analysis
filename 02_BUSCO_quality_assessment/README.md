# Bipolaris Genome and Protein BUSCO Quality Assessment

A reproducible workflow to prepare accession-matched genome/protein FASTA files, install BUSCO 6.1.0, download the `pleosporales_odb12` lineage, run genome and protein BUSCO analyses sequentially, retain successful TXT/JSON summaries, preserve failed outputs/logs, and create combined Excel/TSV reports.

## 1. Expected source layout

The preparation script expects accession-specific folders:

```text
Bipolaris/
├── GCA_000338975.1/
│   ├── *_genomic.fna
│   └── proteins.faa
├── GCA_000338995.1/
│   ├── *_genomic.fna
│   └── proteins.faa
└── BUSCO/
```

It searches recursively outside the `BUSCO` folder. For each accession folder, it prefers one FASTA whose name contains `genomic`. If none exists, it accepts a single non-CDS/non-RNA/non-protein FASTA candidate. Multiple candidates are reported for manual review. Protein input is a file named `proteins.faa`.

Prepared inputs:

```text
BUSCO/
├── Genome input/
│   ├── GCA_000338975.1.fasta
│   └── GCA_000338995.1.fasta
└── Protein input/
    ├── GCA_000338975.1.faa
    └── GCA_000338995.1.faa
```

Genome and protein files must use the same accession basename to pair correctly.

## 2. Install BUSCO 6.1.0

In WSL Ubuntu, initialize Conda (adjust if Miniconda is elsewhere):

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda env create -f environment/busco_environment.yml
conda activate busco_latest
busco --version
which busco
```

Expected: `BUSCO 6.1.0`. If the environment already exists, just activate it.

## 3. Project path

Default project path used in scripts:

```text
/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO
```

Windows equivalent:

```text
D:\\Bipolaris_genomes\\Bipolaris\\BUSCO
```

If your project is elsewhere, edit `BASE` in the shell scripts and Python report script. Copy this repository's `scripts` folder into the project's `BUSCO/scripts` folder, or run the scripts from this repository; the output location is set by `BASE`.

Create project folders:

```bash
BASE="/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO"
mkdir -p "$BASE/Genome input" "$BASE/Protein input" \
 "$BASE/BUSCO_RESULTS/GENOME" "$BASE/BUSCO_RESULTS/PROTEIN" \
 "$BASE/BUSCO_RESULTS/FAILED" "$BASE/BUSCO_RESULTS/LOGS" \
 "$BASE/BUSCO_RESULTS/DOWNLOADS" "$BASE/results"
```

## 4. Prepare input files

Edit `SOURCE_DIR` near the top of `scripts/01_prepare_inputs.py` if necessary. The default is the Windows path `D:\\Bipolaris_genomes\\Bipolaris`.

Run this preparation script using Windows Python:

```bash
python scripts/01_prepare_inputs.py
```

Or from WSL, change `SOURCE_DIR` in the script to `/mnt/d/Bipolaris_genomes/Bipolaris` and run:

```bash
python scripts/01_prepare_inputs.py
```

It copies one genome and one `proteins.faa` per accession, reports matching accessions, missing files and multiple genome candidates, and writes `input_preparation_report.tsv`. Review any multiple-candidate entries. If your 74 `.fasta` and 74 `.faa` files are already prepared, skip this step.

## 5. Download the lineage

```bash
conda activate busco_latest
bash scripts/02_download_lineage.sh
```

Only `pleosporales_odb12` is downloaded. Default location: `~/busco_downloads/lineages/pleosporales_odb12`.

## 6. Run the complete workflow

```bash
conda activate busco_latest
bash scripts/06_run_complete_pipeline.sh
```

The master script verifies BUSCO and the lineage, checks input files, runs all genome files first, then all protein files, and generates the combined reports. Default is 8 threads per sample. Change `THREADS=8` in the two run scripts to change CPU allocation.

Individual stages:

```bash
bash scripts/03_run_genomes.sh
bash scripts/04_run_proteins.sh
python scripts/05_summarize_results.py
```

## 7. BUSCO commands

Genome mode:

```bash
busco -i sample.fasta -m genome -l pleosporales_odb12 -o sample_name -c 8
```

Protein mode:

```bash
busco -i sample.faa -m proteins -l pleosporales_odb12 -o sample_name -c 8
```

## 8. Output structure

```text
BUSCO/
├── Genome input/
├── Protein input/
├── BUSCO_RESULTS/
│   ├── GENOME/ACCESSION/short_summary*.txt + short_summary*.json
│   ├── PROTEIN/ACCESSION/short_summary*.txt + short_summary*.json
│   ├── FAILED/          # failed/incomplete full outputs
│   ├── LOGS/            # one log per sample
│   └── DOWNLOADS/
└── results/
    ├── Bipolaris_BUSCO_Genome_Protein_Summary.xlsx
    └── Bipolaris_BUSCO_Genome_Protein_Summary.tsv
```

For successful runs, detailed output files are permanently deleted after the TXT and JSON summaries have been copied. Failed or incomplete output folders are retained under `FAILED`. Logs are retained. A sample is skipped on rerun only if both non-empty summary files exist.

## 9. BUSCO result interpretation

- C: complete BUSCOs
- S: complete single-copy BUSCOs
- D: complete duplicated BUSCOs
- F: fragmented BUSCOs
- M: missing BUSCOs
- n: total BUSCO groups assessed

The combined report pairs results by accession basename. Project-specific screening labels are:
- `KEEP`: genome and protein C >=95%, and F <=2% for both.
- `REMOVE_FRAGMENTED`: either C <90% or F >5%.
- `REVIEW`: other paired results.
- `INCOMPLETE`: genome or protein result missing.

These are screening thresholds chosen for this project, not universal BUSCO standards. Do not use them as the sole basis for assembly acceptance/rejection.

## 10. Troubleshooting

Conda not found:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate busco_latest
```

Check BUSCO:

```bash
which busco
busco --version
```

Check lineage:

```bash
ls -ld ~/busco_downloads/lineages/pleosporales_odb12
```

Check input counts:

```bash
find "/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO/Genome input" -maxdepth 1 -name '*.fasta' | wc -l
find "/mnt/d/Bipolaris_genomes/Bipolaris/BUSCO/Protein input" -maxdepth 1 -name '*.faa' | wc -l
```

Inspect per-sample logs in `BUSCO_RESULTS/LOGS` and failed outputs in `BUSCO_RESULTS/FAILED`.

## 11. Citation

Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. *Molecular Biology and Evolution*, 38(10), 4647–4654. https://doi.org/10.1093/molbev/msab199

Record the BUSCO version, lineage dataset (`pleosporales_odb12`), and download date in any publication or report.
