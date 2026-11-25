# Prediction output count interpretation

All six AUGUSTUS models were run on genome assembly GCA_000527765.2.

The prediction summary reports gene, transcript, CDS-record and protein-record
counts. These counts agree with the corresponding GFF3 gene/transcript features
and extracted FASTA records.

Important distinction:
- GFF3 `gene` features count predicted genes.
- GFF3 `transcript` features count predicted transcripts.
- GFF3 `CDS` features count individual coding sequence segments.
- A multi-exon transcript can contain multiple CDS features.
- Extracted CDS FASTA records represent complete CDS sequences, not individual
  GFF3 CDS segments.
- Protein FASTA records represent translated transcript products.

Therefore, the total number of GFF3 CDS features must not be reported as the
number of predicted genes.

The prediction summary is:
`prediction_gene_counts.tsv`

These checks establish output consistency and successful prediction execution.
They do not independently establish biological gene-model accuracy.
