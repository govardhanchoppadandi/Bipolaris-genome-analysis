# Bipolaris AUGUSTUS gene-prediction models and workflows

This repository contains six custom AUGUSTUS species parameter sets, prediction workflows, accession manifests, training-provenance records, prediction outputs, and comparative quality-control results for *Bipolaris* genomes.

## Available models

| Version | AUGUSTUS species name |
|---|---|
| V1 | `bipolaris_nd90pr` |
| V2 | `bipolaris_nd90pr_v2` |
| V3 | `bipolaris_multi_v3_all21` |
| V4 | `bipolaris_multi_v4_advanced` |
| V5 clean | `bipolaris_multi_v5_clean` |
| V5 deduplicated | `bipolaris_multi_v5_deduplicated` |

See `metadata/model_copy_report.tsv` for the model-file copy audit and `metadata/model_sha256sums.txt` for SHA-256 checksums.

## Predict genes

Activate the environment in which AUGUSTUS and gffread are installed. Set the AUGUSTUS configuration path and run:

```bash
export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"

bash scripts/predict_one_genome.sh \
  bipolaris_nd90pr_v2 \
  /path/to/genome.fna \
  results/accession_name
```

The example path `/path/to/genome.fna` is a placeholder. Replace it with your own genome FASTA path before running the prediction workflow.

## Validation snapshot

All six models were evaluated on the same *Bipolaris victoriae* test assembly, RefSeq accession `GCA_000527765.2`. BUSCO 6.1.0 was run in protein mode with the `pleosporales_odb12` lineage dataset (6,066 BUSCO groups).

| Model | Predicted genes | Complete BUSCO (%) | Fragmented (%) | Missing (%) |
|---|---:|---:|---:|---:|
| V1 ND90Pr | 9,994 | 93.7 | 1.1 | 5.2 |
| V2 ND90Pr v2 | 10,196 | 94.9 | 1.3 | 3.8 |
| V3 all21 | 9,931 | 93.8 | 1.1 | 5.1 |
| V4 advanced | 9,684 | 92.5 | 1.2 | 6.3 |
| V5 clean | 9,911 | 93.8 | 1.1 | 5.1 |
| V5 deduplicated | 9,539 | Not available | Not available | Not available |

V2 had the highest complete BUSCO percentage and lowest missing percentage among models with available BUSCO results. It is the provisional primary candidate, not a biologically validated final model.

Additional computational checks include:
- Protein/CDS translation consistency.
- CDS length and reading-frame checks.
- Internal stop-codon investigation and genomic mapping.
- Splice-junction motif and phase consistency audits.
- Gene-structure and protein-length summaries.

See `results/VALIDATION_REPORT.md` for methods, detailed results, limitations, and recommended further validation. Supporting QC tables and prediction outputs are included under `test_predictions/`.

## Training provenance

`training_provenance/model_training_status.tsv` distinguishes recovered information from missing provenance. Not every historical script or training accession list has been recovered.

The V3 accession list represents accession folders and does not guarantee 21 independent genomes. The V3 training set included an assembly related to the test assembly, so its BUSCO result is not a fully independent evaluation.

The historical V3 optimization command refers to `bipolaris_multi_v3`, not necessarily the later `bipolaris_multi_v3_all21` model.

## Scope and limitations

These results are computational quality-control evidence. BUSCO completeness, sequence consistency, and canonical splice motifs do not independently establish gene-boundary accuracy or biological correctness.

Further validation should include:
1. BUSCO analysis of V5 deduplicated.
2. Functional-domain annotation of all six protein sets.
3. Independent RNA-seq splice-junction support.
4. Comparison against trusted reference annotations on independent assemblies.
5. Ortholog conservation and domain-architecture analysis.
6. Recovery or clearer documentation of missing training provenance.
7. Testing on additional independent *Bipolaris* genomes.

The test assembly is not redistributed in this repository; obtain it from its source accession.

## Citation

Please cite the archived, versioned repository when available, together with AUGUSTUS, BUSCO, and the source genome or annotation accessions. Add author information and a repository DOI when the release is ready.

## License and reuse

The project-authored source code in `scripts/` and
`original_scripts/` is licensed under the GNU General Public
License v3.0 only (GPL-3.0-only).

The complete license is provided in
`LICENSES/GPL-3.0-only.txt`. See `LICENSES/README.md` for
the scope of this license.

The trained AUGUSTUS parameter files in `models/`, prediction
outputs, reference-derived data, and other third-party materials
are not automatically covered by this code license. Their
redistribution and reuse remain subject to applicable ownership,
attribution, and licensing requirements.

Please cite this repository and the relevant reference assemblies
and software when using the models or workflows.
