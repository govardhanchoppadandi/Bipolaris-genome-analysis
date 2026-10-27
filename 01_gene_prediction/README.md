# Bipolaris Gene Prediction Using AUGUSTUS

## 1. Project overview

This directory contains the initial gene-prediction workflow developed for the **Bipolaris Genome Analysis** project.

The objective is to predict gene structures and derive predicted protein and coding DNA sequences from genome assemblies belonging to *Bipolaris* species.

The workflow uses AUGUSTUS, a gene-prediction program, with a custom species model trained using the annotated ND90Pr reference genome. The program `gffread` is used to extract predicted protein and coding DNA sequences from the predicted gene annotations.

The pipeline is designed for Ubuntu Linux, including Ubuntu 24.04 under Windows Subsystem for Linux 2 (WSL2), with software installed through Conda.

**Project status:** Initial workflow under development. Model optimization, validation, and prediction quality must be reviewed before the results are used for biological interpretation.

---

## 2. Workflow objectives

The workflow covers the following stages:

1. Obtain an annotated reference genome for model training.
2. Prepare training examples from the reference genome annotation.
3. Split the training examples into training and test sets.
4. Train a custom AUGUSTUS species model.
5. Optimize model parameters.
6. Validate the trained model using held-out test data.
7. Predict genes in additional *Bipolaris* genome assemblies.
8. Extract predicted protein and coding DNA sequences.
9. Generate logs and reports documenting the batch prediction process.

---

## 3. Reference genome used for training

The reference selected for this project is *Bipolaris sorokiniana* ND90Pr, also associated with the assembly name Cocsa1 and the name *Cochliobolus sativus* in some historical resources.

| Attribute | Information |
|---|---|
| Organism | *Bipolaris sorokiniana* |
| Strain | ND90Pr |
| Assembly name | Cocsa1 |
| NCBI RefSeq assembly accession | GCF_000338995.1 |
| Purpose | Training the custom AUGUSTUS model |
| AUGUSTUS model name | `bipolaris_nd90pr` |

Reference resources:

- NCBI Assembly: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000338995.1/
- Ensembl Fungi reference information: https://jun2026-fungi.ensembl.org/Bipolaris_sorokiniana_nd90pr_gca_000338995/Info/Index

The ND90Pr genome is the **training reference**. It is not substituted for the genome assemblies being analyzed. During batch prediction, each accession's own genome FASTA is supplied to AUGUSTUS.

The custom model is intended as a starting point for related *Bipolaris* genomes. Prediction performance may differ between species and assemblies.

---

## 4. Reference files to download

Download the reference genome sequence and its corresponding gene annotation from NCBI or another appropriate authoritative source.

The training workflow requires:

| File | Purpose |
|---|---|
| `GCF_000338995.1_Cocsa1_genomic.fna` | Reference genome DNA sequence in FASTA format |
| `genomic.gff` | Reference gene annotation in GFF format |

The following files may also be useful for independent comparison or downstream analyses, but they are not required by the included training script:

| File | Purpose |
|---|---|
| `protein.faa` | Reference proteins |
| `cds_from_genomic.fna` | Reference coding DNA sequences |

The filenames supplied by a download service may differ. If the names differ, use the actual filenames in the commands or rename the files consistently.

### Download procedure

1. Open the NCBI Assembly link above.
2. Locate the download options for assembly `GCF_000338995.1`.
3. Download the genomic sequence in FASTA format.
4. Download the corresponding gene annotation in GFF or GFF3 format.
5. Keep the genome FASTA and annotation from the same assembly and annotation release.
6. Extract the downloaded files if they are compressed.
7. Verify that both files exist and are non-empty before training.

Example reference directory:

```text
ND90Prt/
└── GCF_000338995.1/
    ├── GCF_000338995.1_Cocsa1_genomic.fna
    ├── genomic.gff
    ├── protein.faa
    └── cds_from_genomic.fna
```

The protein and CDS files are shown as optional reference files.

**Important:** Do not use a GFF annotation from one assembly with a genome FASTA from another assembly. The annotation sequence identifiers and coordinates must correspond to the genome sequence used for training.

---

## 5. Input genome directory structure

The prediction workflow processes accession directories located directly inside each input root.

The current project uses these Windows directories:

```text
D:\Bipolaris_genomes\
├── ARI genomes\
│   ├── GCA_accession_1\
│   │   └── genome_genomic.fna
│   ├── GCA_accession_2\
│   │   └── genome_genomic.fna
│   └── ...
│
└── Bipolaris\
    ├── GCA_accession_3\
    │   └── genome_genomic.fna
    ├── GCA_accession_4\
    │   └── genome_genomic.fna
    └── ...
```

Under WSL2, these directories are accessed as:

```text
/mnt/d/Bipolaris_genomes/ARI genomes/
/mnt/d/Bipolaris_genomes/Bipolaris/
```

At the time of initial setup, these roots contained 9 accession folders under `ARI genomes` and 65 under `Bipolaris`.

These numbers describe the initial directory inventory and may change as genomes are added or removed.

### Input requirements

Each accession directory should contain a genome assembly FASTA file directly inside it. The pipeline searches for files ending in `.fna`, `.FNA`, `.fasta`, `.FASTA`, `.fa`, or `.FA`.

Nested assemblies and FASTA files in deeper subdirectories are not automatically searched by the batch script.

---

## 6. Repository structure

The relevant project files are organized as follows:

```text
Bipolaris-genome-analysis/
│
├── README.md
│
└── 01_gene_prediction/
    ├── README.md
    ├── environment.yml
    ├── LICENSE
    ├── CITATION.cff
    ├── .gitignore
    │
    ├── scripts/
    │   ├── train_bipolaris_model.sh
    │   ├── validate_bipolaris_model.sh
    │   └── run_augustus_all_bipolaris.sh
    │
    └── docs/
        └── OUTPUTS.md
```

The scripts perform the following tasks:

- `train_bipolaris_model.sh`: prepares training data, trains the custom model, and attempts parameter optimization.
- `validate_bipolaris_model.sh`: predicts on the held-out test set and saves validation output.
- `run_augustus_all_bipolaris.sh`: runs gene prediction and sequence extraction for accession folders in both input roots.
- `environment.yml`: specifies the main Conda packages.
- `docs/OUTPUTS.md`: describes generated output files.

Genome assemblies, training intermediates, local model files, and large prediction results are not included in this GitHub source-code release.

---

## 7. Software requirements

The workflow requires:

- Conda or Mamba
- AUGUSTUS
- `gffread`
- Perl
- AUGUSTUS helper scripts, including `gff2gbSmallDNA.pl`, `randomSplit.pl`, `new_species.pl`, `etraining`, and `optimize_augustus.pl`

The included `environment.yml` specifies the main software dependencies. After creating or activating the environment, verify that the AUGUSTUS helper scripts are available in the installation.

### Create the environment

From the repository directory, run:

```bash
conda env create -f environment.yml
conda activate bipolaris
```

If the `bipolaris` environment already exists, activate it instead of recreating it.

For a new Ubuntu terminal where Conda is not initialized:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate bipolaris
```

Configure the AUGUSTUS parameter directory:

```bash
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"
```

Verify the main programs:

```bash
which augustus
which gffread
which gff2gbSmallDNA.pl
which randomSplit.pl
which new_species.pl
which etraining
which optimize_augustus.pl
```

All required programs must be available before the training workflow is run.

---

## 8. Training the custom AUGUSTUS model

The custom model is named:

```text
bipolaris_nd90pr
```

The training workflow converts the reference annotation and genome into training examples, divides the examples into training and test sets, and uses AUGUSTUS training tools to fit model parameters.

### Run the training script

Activate the Conda environment and run:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate bipolaris

export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

bash 01_gene_prediction/scripts/train_bipolaris_model.sh \
  "/mnt/d/Bipolaris_genomes/ND90Prt/GCF_000338995.1" \
  "/mnt/d/Bipolaris_genomes/ND90Pr_training" \
  "GCF_000338995.1_Cocsa1_genomic.fna" \
  "genomic.gff"
```

This command assumes that the training script is run from the root of the repository and that the reference files have the names shown above.

If the script is run from inside `01_gene_prediction`, adjust the script path accordingly.

### Training procedure

The training script is intended to:

1. Check that the reference FASTA and GFF annotation are available.
2. Convert annotated genes into GenBank training format using `gff2gbSmallDNA.pl`.
3. Split the examples using `randomSplit.pl`.
4. Create the species model if it does not already exist.
5. Fit model parameters using `etraining`.
6. Run `optimize_augustus.pl` to optimize the model parameters.

Optimization may take considerable time. Run only one optimization process for a given model at a time.

Logs are written to the training directory so that warnings and errors can be reviewed.

### Expected training directory

```text
ND90Pr_training/
├── ND90Pr.genes.gb
├── ND90Pr.genes.gb.train
├── ND90Pr.genes.gb.test
├── gff2gbSmallDNA.log
├── etraining.log
├── optimize.log
├── ND90Pr.validation.out
└── ND90Pr.validation.log
```

The validation files appear after validation has been run. Some files may not exist if an earlier stage failed.

The trained model parameters are stored in the AUGUSTUS configuration directory, normally:

```text
$CONDA_PREFIX/config/species/bipolaris_nd90pr/
```

The exact contents of that directory depend on the installed AUGUSTUS version and the training steps completed.

---

## 9. Model validation

After optimization has finished, validate the model on the held-out test set.

Run:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate bipolaris

export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

bash 01_gene_prediction/scripts/validate_bipolaris_model.sh \
  "/mnt/d/Bipolaris_genomes/ND90Pr_training"
```

The validation script writes:

```text
ND90Pr.validation.out
ND90Pr.validation.log
```

Review the exit status, the generated prediction output, and the log warnings.

A successful command indicates that AUGUSTUS executed successfully; it does not by itself establish the accuracy of every predicted gene.

The held-out test set provides a basic check of model behavior on reference-derived examples. It does not guarantee equal accuracy on all other *Bipolaris* species.

---

## 10. Batch prediction across genome accessions

After model optimization and validation have been reviewed, run the batch prediction script.

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate bipolaris

export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"
export AUGUSTUS_SPECIES="bipolaris_nd90pr"

bash 01_gene_prediction/scripts/run_augustus_all_bipolaris.sh \
  "/mnt/d/Bipolaris_genomes/ARI genomes" \
  "/mnt/d/Bipolaris_genomes/Bipolaris"
```

The script processes accession directories in both input roots.

For each accession, the workflow:

1. Checks whether a complete set of non-empty outputs already exists.
2. Searches for a genome FASTA file directly inside the accession folder.
3. Runs AUGUSTUS using the `bipolaris_nd90pr` model.
4. Writes a GFF3 file with predicted gene structures.
5. Runs `gffread` to extract protein and CDS sequences.
6. Records the result in a tab-separated report.
7. Continues to the next accession if a particular accession fails.

### Skip behavior

An accession is skipped only when all three expected output files exist and are non-empty:

```text
genes.gff3
proteins.faa
cds.fna
```

The script does not independently determine which model was used to generate existing files. If previous outputs were generated with a different model and need to be replaced, back them up and remove the old outputs before rerunning.

---

## 11. Output files

For each successfully processed genome, the script creates:

| File | Description |
|---|---|
| `genes.gff3` | Gene structures predicted by AUGUSTUS |
| `proteins.faa` | Predicted amino-acid sequences extracted by `gffread` |
| `cds.fna` | Predicted coding DNA sequences extracted by `gffread` |

Each input root also receives:

| File | Description |
|---|---|
| `gene_prediction_report.tsv` | Status of each accession processed |
| `augustus_prediction.log` | Combined error and warning output from prediction and extraction |

Possible report statuses include:

- `COMPLETED`
- `SKIPPED`
- `FAILED`
- `NO_FASTA`

A `COMPLETED` status indicates that the prediction and sequence extraction commands completed and produced non-empty output files. It does not establish biological correctness.

See `docs/OUTPUTS.md` for further information.

---

## 12. Monitoring a running analysis

The batch script prints the accession currently being processed.

To inspect the log for the ARI genomes:

```bash
tail -f "/mnt/d/Bipolaris_genomes/ARI genomes/augustus_prediction.log"
```

To inspect the current report:

```bash
cat "/mnt/d/Bipolaris_genomes/ARI genomes/gene_prediction_report.tsv"
```

To inspect the report for the other input root:

```bash
cat "/mnt/d/Bipolaris_genomes/Bipolaris/gene_prediction_report.tsv"
```

To check whether optimization is still running:

```bash
pgrep -af optimize_augustus.pl
```

Wait for optimization to finish before treating the model as fully optimized.

---

## 13. Scientific considerations and limitations

- The model is trained using the ND90Pr / Cocsa1 reference annotation.
- The same model is applied to multiple genome assemblies, so prediction performance may vary with species, evolutionary distance, assembly quality, and annotation conventions.
- The workflow does not incorporate RNA-seq or transcript evidence.
- The included pipeline does not perform functional annotation of predicted proteins.
- Predicted proteins should not be described as experimentally confirmed proteins.
- Gene models and protein sequences should be checked for completeness, plausible coding lengths, and consistency with the genome assembly.
- Further validation using trusted annotations, transcript evidence, or protein homology is recommended where available.
- The training and prediction logs should be retained for reproducibility.

---

## 14. Reproducibility and data management

For reproducible analyses, record:

- Reference assembly accession and annotation source.
- AUGUSTUS version and Conda package versions.
- Custom model name and training parameters.
- Training and validation logs.
- Genome accessions analyzed.
- Batch prediction reports.
- Any manual filtering or post-processing performed.

The repository stores scripts and documentation. Large genome FASTAs, training intermediates, local AUGUSTUS model files, and generated prediction outputs should generally remain outside the Git repository unless there is a specific reason and sufficient repository storage.

The `.gitignore` file helps prevent accidental inclusion of common sequence files and intermediate outputs. Review `git status` before committing.

---

## 15. Planned future extensions

This gene-prediction stage provides a foundation for additional comparative analyses, which may include:

- Genome assembly quality assessment.
- Comparison of gene counts and predicted proteomes.
- Ortholog identification and gene-family analysis.
- Functional annotation and protein-domain analysis.
- Phylogenomic analysis.
- Comparative analysis of candidate pathogenicity-associated genes.

These are planned extensions and should be added only when the corresponding workflows have been implemented and evaluated.

---

## 16. Citation and acknowledgement

Please cite the AUGUSTUS software and the source of the reference assembly when publishing results generated with this workflow.

Reference assembly information:

NCBI RefSeq assembly: `GCF_000338995.1` (Cocsa1; ND90Pr).

The scripts and documentation in this directory are part of the Bipolaris Genome Analysis project.

