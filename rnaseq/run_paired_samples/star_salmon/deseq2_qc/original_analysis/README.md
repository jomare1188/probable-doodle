# Original strict analysis (superseded)

Scripts and results of the first DEA / GO / KEGG analysis of this run (147 up / 82 down DEGs),
kept for reference. They were replaced in `../` because of two problems:

1. `deseq2.dds.RData` (and the nf-core PCA / sample-distance files next to it) came from an
   earlier run: 30,491 rows with `gene-` / `rna-` IDs instead of the 14,137 `DIATSA_LOCUS` genes
   of this run's quantification. The RUVs k = 1 factor used in the DEA was estimated from it.
   The files in `../` were regenerated with nf-core/rnaseq 3.21.0 `bin/deseq2_qc.r` on
   `../../salmon.merged.gene_counts_length_scaled.tsv`, exactly as the pipeline does.
2. `ruv.r` filtered topGO p-values as text (`filter(Classic < 0.05)` on `GenTable` output), which
   dropped every term printed in scientific notation, i.e. the most significant ones
   (e.g. *defense response to bacterium*, p = 3e-10, in the down-regulated genes).

The corrected analysis keeps all 147 original up-regulated genes and 80 of the 82 down-regulated
genes (266 up / 147 down in total). The KEGG results here were also computed with an older KEGG
release than the current ones.
