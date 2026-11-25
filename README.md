# Comprehensive Genome Analysis and Comparative Genomics of *Bipolaris* Species

This repository provides a comprehensive collection of bioinformatics workflows, computational scripts, trained gene-prediction models, supporting datasets, analysis outputs, and documentation for genome analysis and comparative genomics of *Bipolaris* species.

The repository brings together genome annotation, gene prediction, genome quality assessment, functional characterization, phylogenetic analysis, comparative genomics, gene-family evolution, repeat analysis, and other genomic investigations within a structured and reproducible computational framework.

All major analyses are organized into dedicated modules. Each module contains the relevant computational resources, including scripts, software environments, input information, model files where applicable, results, quality-control reports, and detailed explanations of the analytical procedures.

## Study Overview

The objective of this project is to establish an integrated genomic analysis framework for investigating the genome organization, gene content, evolutionary relationships, functional characteristics, and genomic diversity of *Bipolaris* species.

The project incorporates multiple complementary computational approaches, beginning with genome annotation and completeness assessment and extending to comparative genomics, evolutionary analysis, repetitive DNA characterization, and functional interpretation.

Particular emphasis is placed on custom AUGUSTUS gene-prediction model development, systematic evaluation of predicted gene sets, and the organization of reproducible analysis workflows.

## Overall Computational Workflow

The computational framework covers the following major areas:

1. Custom gene prediction using AUGUSTUS.
2. Genome assembly and predicted-protein completeness assessment using BUSCO.
3. Genome assembly quality assessment.
4. Functional annotation of predicted proteins.
5. Phylogenetic analysis of *Bipolaris* and related fungal species.
6. Orthology identification and gene-family analysis.
7. Gene-family expansion and contraction analysis using CAFE5.
8. Non-coding RNA identification and characterization.
9. Transposable-element identification and analysis.
10. Repeat-induced point mutation (RIP) analysis.
11. Transposable-element insertion analysis using McClintock.
12. KEGG orthology and pathway analysis.
13. Genome-wide visualization and Circos analysis.
14. Telomere identification and characterization.
15. Candidate effector protein analysis.

The modules are designed to maintain the connection between input data, computational procedures, generated results, and subsequent biological interpretation.

## Repository Structure

```text
Bipolaris-genome-analysis/
├── 01_AUGUSTUS_Gene_Prediction/
├── 02_BUSCO_quality_assessment/
└── README.md
```

Each numbered directory represents a distinct analysis component. Supporting data, selected results, and supplementary materials are organized separately where appropriate.

## Analysis Modules

| Module | Analysis | Main scope |
|---|---|---|
| 01 | AUGUSTUS Gene Prediction | Custom gene-prediction models, training, prediction, model comparison, and quality control |
| 02 | BUSCO Quality Assessment | Completeness assessment of genome assemblies and predicted proteins |
| 03 | Genome Assembly Quality | Genome assembly statistics and quality evaluation |
| 04 | Functional Annotation | Functional characterization of predicted fungal proteins |
| 05 | Phylogeny | Phylogenetic reconstruction and evolutionary relationships |
| 06 | Orthology and Gene Family Analysis | Ortholog identification and comparative gene-family analysis |
| 07 | CAFE5 Gene Family Evolution | Gene-family expansion and contraction |
| 08 | ncRNA | Non-coding RNA identification and characterization |
| 09 | Transposable Elements | Repeat and transposable-element characterization |
| 10 | RIP Analysis | Repeat-induced point mutation analysis |
| 11 | McClintock | Transposable-element insertion analysis |
| 12 | KEGG | KEGG orthology and pathway analysis |
| 13 | Genome Circos | Genome-wide circular visualization |
| 14 | Telomere Analysis | Identification and characterization of telomeric regions |
| 15 | Effector Prediction | Candidate fungal effector protein analysis |

---

## 01. AUGUSTUS Gene Prediction

The AUGUSTUS module contains custom gene-prediction models developed for *Bipolaris* genome analysis.

Six model versions are included:

- `bipolaris_nd90pr`
- `bipolaris_nd90pr_v2`
- `bipolaris_multi_v3_all21`
- `bipolaris_multi_v4_advanced`
- `bipolaris_multi_v5_clean`
- `bipolaris_multi_v5_deduplicated`

The module maintains the model files and associated computational resources required to understand, reproduce, and evaluate the gene-prediction procedures.

Its contents include:

- Trained AUGUSTUS species-model files.
- Training accession information.
- Available training datasets and provenance records.
- Original and subsequent computational scripts.
- Gene-prediction commands and supporting procedures.
- GFF3 gene-structure annotations.
- Predicted coding DNA sequences (CDS).
- Predicted protein sequences.
- Test prediction results.
- Model metadata and integrity records.
- Sequence-quality assessment reports.
- Comparative model-evaluation results.
- Relevant software information, references, and licensing documentation.

### Model development and evaluation

The six models represent different stages and approaches in custom gene-prediction model development.

Their evaluation includes comparisons of predicted gene numbers, transcript numbers, CDS and protein records, BUSCO completeness, sequence integrity, translation consistency, and other available quality-control measures.

The module also documents the available training provenance and distinguishes verified training information from information that could not be recovered.

Detailed descriptions of each model, its development history, scripts, parameters, prediction outputs, evaluation results, and limitations are provided in the dedicated AUGUSTUS documentation.

**Module documentation:** `01_AUGUSTUS_Gene_Prediction/README.md`

## 02. BUSCO Quality Assessment

The BUSCO module provides computational procedures for assessing the completeness of genome assemblies and predicted protein sets.

The workflow includes:

- Input preparation.
- BUSCO lineage dataset acquisition.
- Genome-based completeness assessment.
- Protein-based completeness assessment.
- Result collection and summarization.
- Automated execution of the available analysis pipeline.
- Software environment configuration.

BUSCO results provide information about complete single-copy, complete duplicated, fragmented, and missing conserved orthologs.

The module supports the evaluation of genome completeness and the comparison of predicted gene sets generated using different gene-prediction models.

**Software environment:** BUSCO 6.1.0 and Python 3.11, as specified in the existing environment configuration.

**Module documentation:** `02_BUSCO_quality_assessment/README.md`

## 03. Genome Assembly Quality

This module organizes the computational procedures used to assess genome assembly characteristics.

The analysis framework covers assembly statistics, sequence-length distributions, contiguity, genome size, and other relevant assembly-quality indicators.

Input information, commands, software requirements, generated reports, and interpretation of assembly metrics are documented within the module.

## 04. Functional Annotation

This module focuses on functional characterization of predicted protein sequences.

The workflow is intended to integrate sequence similarity searches, protein-domain identification, functional annotation, Gene Ontology assignments, and related annotation resources where applicable.

The module documents database requirements, software versions, annotation procedures, and generated functional annotation tables.

## 05. Phylogeny

The phylogenetic analysis module organizes the computational procedures used to investigate evolutionary relationships among *Bipolaris* species and related fungal taxa.

The documentation covers sequence preparation, ortholog or marker selection, alignment, phylogenetic inference, tree processing, visualization, and interpretation.

## 06. Orthology and Gene Family Analysis

This module focuses on comparative analysis of predicted protein sets and identification of orthologous groups.

It documents genome and protein dataset preparation, orthology inference, gene-family classification, comparative analysis, and supporting outputs.

## 07. CAFE5 Gene Family Evolution

The CAFE5 module addresses the evolutionary dynamics of gene families across the selected fungal species.

The analysis framework includes gene-family count preparation, species-tree requirements, configuration, evolutionary modelling, and identification of gene-family expansion and contraction patterns.

## 08. Non-coding RNA Analysis

This module covers the identification and characterization of non-coding RNA sequences in fungal genomes.

The relevant workflows, databases, software requirements, input sequences, and annotation outputs are maintained within the module.

## 09. Transposable Element Analysis

This module focuses on the identification and characterization of transposable elements and repetitive genomic regions.

The analysis includes repeat identification, classification, genomic distribution, and summary of transposable-element characteristics.

## 10. Repeat-Induced Point Mutation Analysis

This module investigates repeat-induced point mutation (RIP), a genome defence-associated process that can influence the evolution and composition of repetitive DNA in fungi.

The analysis documentation covers sequence preparation, RIP-related computational procedures, identification of relevant patterns, and interpretation of the results.

## 11. McClintock Transposable Element Analysis

The McClintock module organizes procedures for comparative transposable-element insertion analysis.

It documents the input requirements, supported tools, execution procedures, and interpretation of insertion-related outputs.

## 12. KEGG Pathway Analysis

This module focuses on KEGG orthology assignment and pathway-level interpretation of predicted fungal proteins.

The documentation covers protein identifiers, KO assignments, annotation resources, pathway mapping, and the organization of resulting tables and reports.

## 13. Genome Circos Visualization

This module contains procedures for visualizing genome-wide characteristics and selected genomic features using circular genome representations.

The workflow documents the preparation of genomic tracks, input tables, visualization parameters, and final graphical outputs.

## 14. Telomere
