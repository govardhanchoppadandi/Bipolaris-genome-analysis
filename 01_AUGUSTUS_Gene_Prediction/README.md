# Custom AUGUSTUS Gene Prediction Models and Evaluation for *Bipolaris* Species

## Project Overview

This module contains six custom AUGUSTUS species parameter sets developed for gene prediction in *Bipolaris* genomes, together with the associated training scripts, prediction workflows, accession information, model metadata, test predictions, and comparative quality-control analyses.

The project focuses on developing and evaluating custom gene-prediction models using different training strategies and datasets. The six models are maintained separately to preserve their individual parameter sets, development histories, and prediction outputs.

The repository also provides computational procedures for examining predicted gene structures, coding sequences (CDS), protein sequences, translation consistency, internal stop codons, splice-junction characteristics, and BUSCO completeness.

The objective is to provide a transparent and reproducible resource for fungal gene prediction, model comparison, and subsequent comparative-genomic analyses.

## Objectives

The major objectives of this module are:

1. Develop custom AUGUSTUS species models for *Bipolaris* genomes.
2. Preserve multiple model versions and their associated parameter files.
3. Maintain available training accession information and model-development provenance.
4. Provide original and subsequent training and prediction scripts.
5. Generate gene-structure predictions in GFF3 format.
6. Extract predicted CDS and protein sequences.
7. Compare the six models using a common test genome.
8. Evaluate predicted protein completeness using BUSCO.
9. Examine CDS and protein sequence consistency.
10. Investigate internal stop codons and their genomic locations.
11. Examine CDS reading-frame characteristics and splice-junction features.
12. Maintain model-integrity records and SHA-256 checksums.
13. Identify limitations and requirements for further independent validation.

## Available AUGUSTUS Models

Six custom AUGUSTUS models are included.

| Version | Model name | Description |
|---|---|---|
| V1 | `bipolaris_nd90pr` | Initial ND90Pr-based model |
| V2 | `bipolaris_nd90pr_v2` | Revised ND90Pr-based model |
| V3 | `bipolaris_multi_v3_all21` | Multi-accession model |
| V4 | `bipolaris_multi_v4_advanced` | Advanced multi-genome model |
| V5 | `bipolaris_multi_v5_clean` | Cleaned training-dataset model |
| V5 | `bipolaris_multi_v5_deduplicated` | Deduplicated training-dataset model |

Each model is stored in its own directory under `models/`.

The model directories contain AUGUSTUS parameter and probability files, including exon, intron, intergenic, metaparameter, and weight-matrix information.

The original parameter files and associated model files are retained to preserve the available model versions.

## Repository Structure

```text
01_AUGUSTUS_Gene_Prediction/
│
├── LICENSES/
│   ├── GPL-3.0-only.txt
│   └── README.md
│
├── accessions/
│   ├── V1_ND90Pr_reference.txt
│   ├── V3_all21_accessions.txt
│   └── accession_provenance.tsv
│
├── metadata/
│   ├── model_copy_report.tsv
│   ├── model_sha256sums.txt
│   ├── original_script_copy_report.tsv
│   └── software_environment.txt
│
├── models/
│   ├── bipolaris_nd90pr/
│   ├── bipolaris_nd90pr_v2/
│   ├── bipolaris_multi_v3_all21/
│   ├── bipolaris_multi_v4_advanced/
│   ├── bipolaris_multi_v5_clean/
│   └── bipolaris_multi_v5_deduplicated/
│
├── original_scripts/
│   └── Original training and prediction scripts
│
├── results/
│   ├── model_validation/
│   ├── prediction_test/
│   ├── sequence_quality/
│   ├── protein_BUSCO_comparison.tsv
│   └── VALIDATION_REPORT.md
│
├── scripts/
│   └── Prediction, validation, and sequence-QC scripts
│
├── test_predictions/
│   ├── cds/
│   ├── proteins/
│   ├── GFF3 prediction files
│   └── Sequence-integrity and gene-structure reports
│
├── training_provenance/
│   └── model_training_status.tsv
│
├── CITATION.cff
├── README.md
├── .gitattributes
└── .gitignore
```

## Model Development and Training Provenance

The six models represent different stages of model development and training.

Available records are maintained in the `accessions/`, `training_provenance/`, and `original_scripts/` directories.

### V1 — ND90Pr

Species model: `bipolaris_nd90pr`

This model is based on the *Bipolaris sorokiniana* ND90Pr reference genome, RefSeq accession `GCF_000338995.1`.

The available documentation describes the ND90Pr-based model-development approach. However, the exact historical execution script and complete command options have not been recovered.

### V2 — ND90Pr version 2

Species model: `bipolaris_nd90pr_v2`

This model represents a revised ND90Pr-based gene-prediction model.

The complete original training accession inventory, split manifest, and historical execution logs have not been recovered. These provenance limitations are retained explicitly rather than inferred.

### V3 — Multi-accession model

Species model: `bipolaris_multi_v3_all21`

This model was developed using a multi-accession training approach.

The available accession information is provided in:

`accessions/V3_all21_accessions.txt`

The accession list represents accession folders and should not automatically be interpreted as 21 independent genomes. Paired assembly accessions may refer to related records.

The V3 training set included an assembly related to the test assembly, `GCA_000527765.2`. Therefore, its evaluation is not a fully independent test.

The recovered historical optimization command is preserved separately in:

`scripts/recovered_historical_v3_optimization_command.sh`

The historical command refers to `bipolaris_multi_v3`, which should not automatically be treated as identical to the final `bipolaris_multi_v3_all21` model.

### V4 — Advanced multi-genome model

Species model: `bipolaris_multi_v4_advanced`

The available development information describes a training set involving 5,250 records and a two-round optimization approach.

The exact complete accession inventory and historical execution logs have not been recovered.

### V5 — Clean model

Species model: `bipolaris_multi_v5_clean`

This model was developed using a cleaned training dataset.

Available development records describe approximately 57,490 GenBank records, with deduplication and test-overlap filtering, followed by two rounds of optimization using eight CPUs.

The complete original source accession inventory and all dataset-construction scripts have not been recovered.

The recovered optimization command is preserved in:

`scripts/recovered_v5_optimization_command.sh`

### V5 — Deduplicated model

Species model: `bipolaris_multi_v5_deduplicated`

This model represents a deduplicated training-data version.

Its parameter files are retained independently. Complete training provenance and a corresponding BUSCO evaluation are not currently available.

### Training provenance records

The file:

`training_provenance/model_training_status.tsv`

summarizes the available model-training information and distinguishes recovered information from missing historical records.

The repository does not claim that all historical training procedures have been completely reconstructed.

## Model Installation and Prediction

### Requirements

The prediction workflow requires:

- AUGUSTUS.
- A configured AUGUSTUS species-parameter directory.
- `gffread` for sequence extraction.
- Bash and the supporting command-line utilities used by the scripts.

The recorded software environment is available in:

`metadata/software_environment.txt`

### AUGUSTUS configuration

Activate the environment containing AUGUSTUS and configure the AUGUSTUS configuration directory.

```bash
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"
```

The custom species models must be available under the configured AUGUSTUS species directory.

The six model directories are provided under `models/`.

### Predicting genes using a custom model

The general prediction script is:

`scripts/predict_one_genome.sh`

Example:

```bash
bash scripts/predict_one_genome.sh \
  bipolaris_nd90pr_v2 \
  /path/to/genome.fna \
  results/accession_name
```

Replace `/path/to/genome.fna` with the actual genome FASTA path.

Replace `accession_name` with the desired output identifier.

The model name must correspond to one of the six installed AUGUSTUS species models.

### Prediction output formats

The test prediction results include:

- GFF3 gene-structure annotations.
- CDS FASTA files.
- Protein FASTA files.

The GFF3 files describe the predicted gene structures and their genomic coordinates.

The CDS files contain predicted coding sequences.

The protein FASTA files contain the corresponding predicted amino acid sequences.

The supplied test predictions are organized under:

`test_predictions/`

## Common Test Genome and Evaluation Design

All six models were evaluated on the same *Bipolaris victoriae* genome assembly.

**Test genome**

- Organism: *Bipolaris victoriae*
- RefSeq assembly accession: `GCA_000527765.2`
- Assembly: `GCA_000527765.2_ANU_Cvic_canu1.6_pilon5`
- Sequence identifier: `AMCY02000001.1`

The same test assembly was used for the six prediction runs to facilitate comparative evaluation.

The test genome itself is not redistributed in this repository. Obtain it from its original accession source.

## Gene Prediction Statistics

The following table summarizes the predicted gene, transcript, CDS, and protein record counts.

| Model | Genes | Transcripts | CDS records | Protein records |
|---|---:|---:|---:|---:|
| V1 ND90Pr | 9,994 | 9,994 | 9,994 | 9,994 |
| V2 ND90Pr v2 | 10,196 | 10,196 | 10,196 | 10,196 |
| V3 all21 | 9,931 | 9,931 | 9,931 | 9,931 |
| V4 advanced | 9,684 | 9,684 | 9,684 | 9,684 |
| V5 clean | 9,911 | 9,911 | 9,911 | 9,911 |
| V5 deduplicated | 9,539 | 9,539 | 9,539 | 9,539 |

The CDS feature count in a GFF3 file can be greater than the gene count because a multi-exon gene may contain multiple CDS features.

The FASTA record counts above represent the corresponding gene/transcript sequence records, not the number of individual GFF3 CDS features.

## BUSCO Evaluation

BUSCO was used to compare the completeness of the predicted protein sets.

**Evaluation configuration**

- BUSCO version: 6.1.0
- Lineage dataset: `pleosporales_odb12`
- Number of BUSCO groups: 6,066
- Analysis mode: Protein
- Test assembly: *Bipolaris victoriae* `GCA_000527765.2`

### BUSCO comparison

| Model | Genes | Complete BUSCO | Single-copy | Duplicated | Fragmented | Missing |
|---|---:|---:|---:|---:|---:|---:|
| V1 ND90Pr | 9,994 | 93.7% | 93.6% | 0.1% | 1.1% | 5.2% |
| V2 ND90Pr v2 | 10,196 | 94.9% | 94.7% | 0.1% | 1.3% | 3.8% |
| V3 all21 | 9,931 | 93.8% | 93.7% | 0.1% | 1.1% | 5.1% |
| V4 advanced | 9,684 | 92.5% | 92.4% | 0.1% | 1.2% | 6.3% |
| V5 clean | 9,911 | 93.8% | 93.7% | 0.1% | 1.1% | 5.1% |
| V5 deduplicated | 9,539 | Not available | Not available | Not available | Not available | Not available |

### Interpretation

Among the models with available BUSCO results, V2 achieved:

- The highest complete BUSCO percentage (94.9%).
- The highest single-copy BUSCO percentage (94.7%).
- The lowest missing BUSCO percentage (3.8%).

V2 is therefore considered the provisional primary candidate based on the currently available completeness results.

This recommendation is provisional. BUSCO completeness alone does not establish the biological accuracy of gene boundaries, exon structures, or all predicted protein sequences.

The V5 deduplicated model requires a corresponding BUSCO evaluation before a complete six-model comparison can be made.

Detailed comparison results are provided in:

`results/protein_BUSCO_comparison.tsv`

## Protein and CDS Sequence Integrity

Additional sequence-level quality-control procedures were conducted to assess the consistency of predicted proteins and CDS sequences.

The analyses include:

- Protein/CDS correspondence.
- Translation consistency.
- CDS length and reading-frame characteristics.
- Internal stop-codon identification.
- Genomic mapping of internal stops.
- Splice-junction motif inspection.
- CDS phase consistency.
- Gene-structure and protein-length summaries.

### Translation consistency

The predicted protein sequences were compared against translations of their corresponding CDS sequences.

AUGUSTUS terminal stop representation was normalized during the comparison.

The corrected translation comparison identified zero mismatches across all protein/CDS pairs in the six model sets.

### CDS integrity

The following table summarizes the sequence-level findings.

| Model | Protein/CDS pairs | Non-triplet CDS records | Internal stop-containing proteins | Translation mismatches |
|---|---:|---:|---:|---:|
| V1 ND90Pr | 9,994 | 0 | 39 | 0 |
| V2 ND90Pr v2 | 10,196 | 0 | 27 | 0 |
| V3 all21 | 9,931 | 2 | 7 | 0 |
| V4 advanced | 9,684 | 1 | 13 | 0 |
| V5 clean | 9,911 | 2 | 12 | 0 |
| V5 deduplicated | 9,539 | 0 | 39 | 0 |

The non-triplet CDS records were retained for review rather than silently removed.

Their details are available in:

`results/sequence_quality/non_triplet_cds_review.tsv`

### Internal stop-codon investigation

A total of 137 internal stop-containing protein records were identified across the six model sets.

The investigation included:

- Identification of internal stop codons.
- Classification of stop-codon types.
- Mapping to genomic coordinates.
- Examination of corresponding CDS features.
- Comparison of recurring genomic loci among models.
- Inspection of splice-junction relationships.

The identified internal stops were canonical stop codons (TAA, TAG, and TGA).

The investigated stop codons were reconstructed from the corresponding genome and GFF3 annotations and were associated with annotated CDS exon junctions.

The corresponding splice-junction motifs were canonical GT-AG in the examined cases.

The CDS phase audit did not identify downstream phase mismatches.

These findings describe the computational annotation and sequence reconstruction results. They do not independently establish the biological correctness of every affected gene model.

### Internal stop recurrence

The internal stop investigation identified 94 unique genomic loci across the model predictions.

Their recurrence was:

| Number of models containing a locus | Number of loci |
|---|---:|
| 5 models | 1 |
| 4 models | 1 |
| 3 models | 10 |
| 2 models | 16 |
| 1 model | 66 |

This comparison helps distinguish recurrent genomic features from model-specific prediction patterns.

### Non-triplet CDS records

Five non-triplet CDS records were identified in V3, V4, and V5 clean.

| Model | Transcript | CDS length |
|---|---|---:|
| V3 all21 | g5734.t1 | 2,141 nt |
| V3 all21 | g9926.t1 | 550 nt |
| V4 advanced | g5586.t1 | 2,192 nt |
| V5 clean | g5713.t1 | 2,141 nt |
| V5 clean | g9907.t1 | 550 nt |

None of these five records was located at a physical contig boundary.

They are retained in the QC report for further review.

## Model File Integrity

The six AUGUSTUS models are preserved as separate model directories.

The model-copy audit records:

- Six model directories.
- Nine model files per directory.
- 54 model files in total.
- Successful copy verification for the recorded files.
- SHA-256 checksum verification.

The relevant records are:

- `metadata/model_copy_report.tsv`
- `metadata/model_sha256sums.txt`

These files provide a record for checking the integrity of the distributed model parameter files.

## Scripts and Computational Resources

The `scripts/` directory contains the computational scripts used for prediction, validation, sequence analysis, and quality control.

The `original_scripts/` directory preserves the available original training, data-preparation, downloading, prediction, and validation scripts.

The scripts are maintained separately so that the original computational records remain distinguishable from later analysis and QC procedures.

The available scripts include:

### Prediction and validation

- `predict_one_genome.sh`
- `validate_prediction_outputs.sh`

### Protein and CDS quality control

- `audit_protein_cds.py`
- `audit_protein_gene_structure.py`
- `audit_internal_stops.py`
- `inspect_translation_mismatches.py`
- `inspect_nontriplet_cds.py`
- `inspect_nontriplet_gff.py`

### Internal stop and splice-junction analysis

- `summarize_internal_stops.py`
- `map_internal_stops_to_genome.py`
- `classify_stop_codon_locations.py`
- `compare_internal_stop_loci.py`
- `inspect_stop_codon_examples.py`
- `audit_stop_codon_phases.py`
- `audit_stop_junction_motifs.py`

### Gene-structure and protein inspection

- `check_gene_boundaries.py`
- `inspect_extreme_proteins.py`

### Recovered historical optimization commands

- `recovered_historical_v3_optimization_command.sh`
- `recovered_v5_optimization_command.sh`

The complete scripts are available in their respective directories. Their detailed usage, arguments, input requirements, and expected outputs should be followed from the individual script contents and supporting documentation.

## Results and Supporting Files

The principal results are maintained in the `results/` and `test_predictions/` directories.

Important reports include:

| File | Description |
|---|---|
| `results/VALIDATION_REPORT.md` | Detailed model validation report |
| `results/protein_BUSCO_comparison.tsv` | Comparative BUSCO results |
| `results/model_validation/prediction_output_integrity.tsv` | Prediction-output integrity information |
| `results/prediction_test/prediction_gene_counts.tsv` | Prediction gene-count summary |
| `results/sequence_quality/non_triplet_cds_review.tsv` | Non-triplet CDS review |
| `test_predictions/sequence_integrity_audit_corrected.tsv` | Corrected sequence-integrity audit |
| `test_predictions/internal_stop_details.tsv` | Internal stop details |
| `test_predictions/internal_stop_genomic_coordinates.tsv` | Genomic coordinates of internal stops |
| `test_predictions/internal_stop_locus_recurrence.tsv` | Recurrence across models |
| `test_predictions/gene_structure_summary.tsv` | Gene-structure summary |
| `test_predictions/protein_length_summary.tsv` | Protein-length summary |
| `test_predictions/stop_codon_phase_audit.tsv` | CDS phase audit |
| `test_predictions/stop_junction_splice_motifs.tsv` | Splice-junction motif results |

The test prediction directory also contains the six GFF3 files and corresponding CDS and protein FASTA files.

## Model Comparison and Current Recommendation

Based on the available prediction counts, BUSCO completeness, and sequence-integrity assessments, V2 (`bipolaris_nd90pr_v2`) is the provisional primary model.

This selection is based on the current computational evidence, particularly its higher complete BUSCO percentage and lower missing BUSCO percentage among the models evaluated.

However, a final biological accuracy assessment requires additional independent evidence.

The six models are retained because different training strategies may produce different gene predictions and may be useful for comparative analyses.

The recommendation should be revisited after additional validation.

## Limitations and Further Validation

The following limitations are recognized:

1. Complete historical training provenance has not been recovered for every model.
2. The V3 training set contains an assembly related to the test assembly, limiting the independence of its evaluation.
3. The BUSCO evaluation for V5 deduplicated is not currently available.
4. Independent RNA-seq splice-junction support has not yet been used as a complete benchmark across all models.
5. Comparison against trusted reference gene annotations on independent assemblies remains necessary.
6. Functional-domain annotation and ortholog conservation should be used as additional evidence.
7. BUSCO completeness and sequence-integrity checks do not independently establish the accuracy of every predicted gene structure.

Further validation should include independent *Bipolaris* genomes, RNA-seq evidence, trusted reference annotations, protein-domain information, and ortholog-based comparisons.

## Citation

Please cite the relevant AUGUSTUS and BUSCO publications, the source genome assembly, and this repository when using the models or computational workflows.

BUSCO reference:

Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. *Molecular Biology and Evolution*, 38(10), 4647–4654. https://doi.org/10.1093/molbev/msab199

Additional software references and version information are maintained in the relevant metadata and environment files.

## License and Reuse

The project-authored source code in:

- `scripts/`
- `original_scripts/`

is licensed under the GNU General Public License v3.0 only (GPL-3.0-only).

The complete license is provided in:

`LICENSES/GPL-3.0-only.txt`

See `LICENSES/README.md` for the scope of the license.

The trained AUGUSTUS parameter files in `models/`, prediction outputs, reference-derived data, and other third-party materials are not automatically covered by this code license. Their redistribution and reuse remain subject to applicable ownership, attribution, and licensing requirements.

Please cite this repository, relevant software, and source genome assemblies when using the models or workflows.

## Author

**Govardhan Choppadandi**

GitHub: https://github.com/govardhanchoppadandi

Repository: https://github.com/govardhanchoppadandi/Bipolaris-genome-analysis
