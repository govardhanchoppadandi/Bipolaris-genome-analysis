# Bipolaris Genome Analysis

## Project overview

This repository contains bioinformatics workflows for genome analysis of *Bipolaris sorokiniana* (syn. *Cochliobolus sativus*).

The project includes genome annotation-related analysis and genome/protein completeness assessment using BUSCO.

## Repository structure

```text
Bipolaris-genome-analysis/
├── 01_gene_prediction/
└── 02_BUSCO_quality_assessment/
    ├── environment/
    ├── results/
    ├── 01_prepare_inputs.py
    ├── 02_download_lineage.sh
    ├── 03_run_genomes.sh
    ├── 04_run_proteins.sh
    ├── 05_summarize_results.py
    ├── 06_run_complete_pipeline.sh
    └── README.md
```

## 1. Gene prediction

See [01_gene_prediction](./01_gene_prediction/) for the gene prediction workflow and related files.

## 2. BUSCO quality assessment

See [02_BUSCO_quality_assessment](./02_BUSCO_quality_assessment/) for the reproducible BUSCO workflow.

The BUSCO analysis uses:

- BUSCO version: 6.1.0
- Lineage dataset: `pleosporales_odb12`
- Genome mode and protein mode
- 74 genome assemblies and 74 corresponding protein datasets

## Reproducibility

Workflow instructions, software environment, scripts, and result organization are documented in the respective analysis folders.

## Organism

*Bipolaris sorokiniana*

## Citation

Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. *Molecular Biology and Evolution*, 38(10), 4647–4654. https://doi.org/10.1093/molbev/msab199
