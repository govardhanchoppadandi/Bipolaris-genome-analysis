# Bipolaris Genome Analysis

A comprehensive and reproducible collection of computational workflows for genome analysis and comparative genomics of *Bipolaris* species.

## Project Overview

This repository organizes gene prediction, genome quality assessment, functional annotation, phylogenetics, orthology, gene-family evolution, repeat analysis, and comparative-genomics workflows into numbered modules.

## Analysis Modules

| Module | Analysis |
|---|---|
| 01 | AUGUSTUS Gene Prediction |
| 02 | BUSCO Quality Assessment |
| 03 | Genome Assembly Quality |
| 04 | Functional Annotation |
| 05 | Phylogeny |
| 06 | Orthology and Gene Family Analysis |
| 07 | CAFE5 Gene Family Evolution |
| 08 | Non-coding RNA Analysis |
| 09 | Transposable Element Analysis |
| 10 | Repeat-Induced Point Mutation (RIP) Analysis |
| 11 | McClintock Transposable Element Analysis |
| 12 | KEGG Pathway Analysis |
| 13 | Genome Circos Visualization |
| 14 | Telomere Analysis |
| 15 | Effector Prediction |

Modules will be populated and documented as their analyses are prepared. Their listing does not imply that every analysis is complete.

## Module 01 — AUGUSTUS Gene Prediction

This module contains six custom AUGUSTUS models:

- `bipolaris_nd90pr`
- `bipolaris_nd90pr_v2`
- `bipolaris_multi_v3_all21`
- `bipolaris_multi_v4_advanced`
- `bipolaris_multi_v5_clean`
- `bipolaris_multi_v5_deduplicated`

It also includes scripts, available training metadata, prediction outputs, quality-control reports, and model-integrity checks.

See [Module 01 README](01_AUGUSTUS_Gene_Prediction/README.md) for detailed instructions.

## Module 02 — BUSCO Quality Assessment

The existing BUSCO module is retained in `02_BUSCO_quality_assessment/`. It documents completeness assessment of genome assemblies and predicted gene sets.

## Reproducibility

Each analysis module should document its software requirements, versions, input data, commands, parameters, outputs, and quality-control procedures.

Large raw datasets may be hosted separately. Refer to the relevant module for accession information and data availability.

## Citation and Licensing

Please cite the relevant software, reference datasets, and repository version when using these workflows or results. Consult each module's licensing documentation because licensing may differ between scripts, models, results, and third-party data.

## Author

Govardhan Choppadandi

GitHub: https://github.com/govardhanchoppadandi
