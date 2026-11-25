# Comparative QC and validation of six Bipolaris AUGUSTUS models

## 1. Scope

Six custom AUGUSTUS species parameter sets were evaluated on the same
Bipolaris victoriae genome assembly, GCA_000527765.2.

The models are:

1. bipolaris_nd90pr
2. bipolaris_nd90pr_v2
3. bipolaris_multi_v3_all21
4. bipolaris_multi_v4_advanced
5. bipolaris_multi_v5_clean
6. bipolaris_multi_v5_deduplicated

This report documents computational execution, output consistency,
BUSCO completeness, sequence integrity, gene structure, and known
provenance limitations.

These analyses provide comparative computational evidence. They do not
by themselves prove biological correctness of every predicted gene.

## 2. Model files and prediction execution

All six model directories were copied successfully. Each contains nine
AUGUSTUS parameter/configuration files: 54 files in total.

The model-copy report is:
`metadata/model_copy_report.tsv`

The model checksums are recorded in:
`metadata/model_sha256sums.txt`

All 54 checksums passed the recorded SHA-256 verification.

All six models completed prediction on GCA_000527765.2 with exit status 0.

| Model | Predicted genes | Transcripts | CDS records | Protein records |
|---|---:|---:|---:|---:|
| bipolaris_nd90pr | 9,994 | 9,994 | 9,994 | 9,994 |
| bipolaris_nd90pr_v2 | 10,196 | 10,196 | 10,196 | 10,196 |
| bipolaris_multi_v3_all21 | 9,931 | 9,931 | 9,931 | 9,931 |
| bipolaris_multi_v4_advanced | 9,684 | 9,684 | 9,684 | 9,684 |
| bipolaris_multi_v5_clean | 9,911 | 9,911 | 9,911 | 9,911 |
| bipolaris_multi_v5_deduplicated | 9,539 | 9,539 | 9,539 | 9,539 |

The GFF3 gene and transcript feature counts agree with the prediction
summary. Extracted protein and complete-CDS FASTA record counts also agree.

A GFF3 CDS feature is an individual coding segment. A multi-exon gene
can have multiple CDS features. Therefore, the total number of GFF3 CDS
features must not be confused with the number of predicted genes.

Detailed counts:
`results/prediction_test/prediction_gene_counts.tsv`

## 3. BUSCO comparison

BUSCO 6.1.0 was run using the pleosporales_odb12 lineage
(6,066 BUSCO groups) on the common test assembly.

| Model | Complete | Single-copy | Duplicated | Fragmented | Missing |
|---|---:|---:|---:|---:|---:|
| bipolaris_nd90pr | 93.7% | 93.6% | 0.1% | 1.1% | 5.2% |
| bipolaris_nd90pr_v2 | 94.9% | 94.7% | 0.1% | 1.3% | 3.8% |
| bipolaris_multi_v3_all21 | 93.8% | 93.7% | 0.1% | 1.1% | 5.1% |
| bipolaris_multi_v4_advanced | 92.5% | 92.4% | 0.1% | 1.2% | 6.3% |
| bipolaris_multi_v5_clean | 93.8% | 93.7% | 0.1% | 1.1% | 5.1% |
| bipolaris_multi_v5_deduplicated | Not available | Not available | Not available | Not available | Not available |

V2 has the highest complete and single-copy BUSCO percentages and the
lowest missing percentage among models with available BUSCO results.
Its fragmented percentage (1.3%) is slightly higher than that of V1,
V3, and V5 clean (1.1%).

BUSCO measures conserved gene-space completeness. It does not directly
measure exact exon boundaries, transcript support, or gene-level
precision and recall.

Source table:
`results/protein_BUSCO_comparison.tsv`

## 4. Protein/CDS sequence integrity

Protein and CDS records were compared for all six models.

| Model | Protein records | Non-triplet CDS | Internal stop codons | Translation mismatches after normalization |
|---|---:|---:|---:|---:|
| bipolaris_nd90pr | 9,994 | 0 | 39 | 0 |
| bipolaris_nd90pr_v2 | 10,196 | 0 | 27 | 0 |
| bipolaris_multi_v3_all21 | 9,931 | 2 | 7 | 0 |
| bipolaris_multi_v4_advanced | 9,684 | 1 | 13 | 0 |
| bipolaris_multi_v5_clean | 9,911 | 2 | 12 | 0 |
| bipolaris_multi_v5_deduplicated | 9,539 | 0 | 39 | 0 |

AUGUSTUS protein output can represent internal stop symbols using
periods. Translation comparisons normalized these symbols before
comparison. No protein/CDS translation mismatches remained after
normalization.

The five non-triplet CDS records are:

- V3: g5734.t1, 2,141 nt; g9926.t1, 550 nt
- V4: g5586.t1, 2,192 nt
- V5 clean: g5713.t1, 2,141 nt; g9907.t1, 550 nt

Coordinate checks showed that none of these five CDS features touches
a physical contig boundary. They are retained as unresolved
non-triplet CDS records requiring review, not labelled as contig-edge
partial genes.

Audit:
`results/sequence_quality/non_triplet_cds_review.tsv`

## 5. Internal stop-codon investigation

Across the six models, 137 internal stop codons were investigated.

The reported findings were:

- All 137 were canonical TAA, TAG, or TGA codons.
- The codons contained no ambiguous bases.
- All were reconstructed from the genomic sequence and GFF3.
- All crossed annotated CDS exon junctions.
- Associated splice junctions had canonical GT-AG motifs.
- The stop-codon phase audit found no downstream phase mismatches.

The loci were grouped into 94 unique genomic locations. Some locations
were shared by multiple models.

These checks support sequence reconstruction and annotation consistency.
Canonical splice motifs and phase consistency do not independently
establish transcript expression or biological validity. RNA-seq or
other transcript evidence is needed.

Relevant audit tables are in `test_predictions/`, including:

- `internal_stop_details.tsv`
- `internal_stop_genomic_coordinates.tsv`
- `internal_stop_locus_recurrence.tsv`
- `stop_codon_phase_audit.tsv`
- `stop_codons_crossing_exon_junctions.tsv`
- `stop_junction_splice_motifs.tsv`
- `translation_mismatch_details.tsv`

## 6. Gene structure and protein length

The gene-structure audit found multi-exon and single-exon predictions
in all six models. Mean exon counts ranged from 2.82 to 3.06 per
transcript.

All six models had zero missing-start records in the reported
gene-structure audit. Missing-stop records were reported for V3
(two), V4 (one), and V5 clean (two); V1, V2, and V5 deduplicated
had zero.

Protein-length summaries identified unusually long predicted proteins
in all six models, including proteins exceeding 8,000 amino acids.
Extreme length alone is not evidence of an incorrect prediction.
These proteins warrant orthology, domain-architecture, and transcript
evidence review.

Relevant tables:

- `test_predictions/gene_structure_summary.tsv`
- `test_predictions/protein_length_summary.tsv`
- `test_predictions/extreme_protein_details.tsv`
- `test_predictions/gene_structure_details.tsv`

## 7. Provisional model assessment

Based on the currently available comparative evidence,
`bipolaris_nd90pr_v2` (V2) is the recommended provisional primary
model.

Reasons:

- Highest complete BUSCO score among models with available BUSCO data.
- Highest single-copy BUSCO score.
- Lowest missing BUSCO percentage.
- Zero non-triplet CDS records in the current sequence-integrity audit.
- Zero normalized protein/CDS translation mismatches.
- Matching GFF3 gene/transcript and extracted FASTA record counts.
- No missing-start or missing-stop records in the reported
  gene-structure audit.

This is a provisional ranking, not proof that V2 has the highest
biological gene-model accuracy.

V3 and V5 clean remain useful alternatives. V5 deduplicated cannot
be fully compared by BUSCO until its BUSCO analysis is available.

## 8. Training provenance and independence limitations

Training provenance is not equally complete for all six models.

- V1: ND90Pr/Cocsa1 reference accession is documented, but the exact
  executed original training script and options were not recovered.
- V2: training accession inventory, split manifest, original script,
  and logs were not recovered.
- V3: an accession list and training script are present. The list
  contains 21 accession records/folders, but paired GCA/GCF records
  may not represent independent genomes. The training set included
  GCA_000527765.1, related to the GCA_000527765.2 test assembly.
  Therefore, the V3 test is not fully independent.
- V4: a 5,250-record sample and two-round optimization are documented,
  but the exact accession list and full execution provenance are
  incomplete.
- V5 clean: 57,490 GenBank training records were reported, with
  deduplication and test-overlap filtering. Full data-construction
  provenance and the complete source accession list are not available.
- V5 deduplicated: alternative model identity is recorded, but its
  complete training provenance is not recovered.

See:

- `training_provenance/model_training_status.tsv`
- `accessions/accession_provenance.tsv`
- `accessions/V1_ND90Pr_reference.txt`
- `accessions/V3_all21_accessions.txt`

Do not describe unavailable scripts, accession lists, or logs as
recovered. The scripts in `original_scripts/` should be interpreted
alongside the provenance tables; their presence does not prove that
each was the exact script executed to create the corresponding final
model.

## 9. Functional-domain analysis

InterProScan/domain annotation has not yet been completed for all six
protein sets in this validation package.

Consequently, no comparative claims about InterPro domain coverage,
Pfam assignments, or domain architecture are made here.

A future analysis should compare the six protein sets using the same
InterProScan version and database release, and report the proportion
of proteins with matches, domain coverage, and conserved domain
architecture. Functional annotation is complementary evidence and
does not by itself validate gene boundaries.

## 10. Further validation required

The following would strengthen biological validation:

1. Run BUSCO for V5 deduplicated under the same settings.
2. Evaluate all six protein sets with the same InterProScan release.
3. Obtain independent RNA-seq splice-junction support.
4. Compare predictions with a trusted reference annotation on an
   independent assembly.
5. Assess ortholog conservation and domain architecture, especially
   for unusually long proteins.
6. Recover or clearly document missing training provenance.
7. Test on additional independent Bipolaris genomes.
8. Where possible, calculate gene-level precision/recall and exon
   boundary agreement against a high-quality reference annotation.

## 11. Overall conclusion

All six custom AUGUSTUS models are available as parameter sets and
completed prediction on the common test assembly. Their prediction
outputs passed the reported record-count consistency checks.
Sequence-integrity, internal-stop, splice-motif, phase, protein-length,
gene-structure, and BUSCO comparisons provide multiple layers of
computational QC.

V2 is the best-supported provisional primary candidate based on the
currently available evidence. The evidence does not establish that
all predicted gene structures are biologically correct. Independent
transcript evidence, additional reference-genome testing, functional
domain analysis, and improved training provenance remain important
next steps.
