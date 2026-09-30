# Diatraea control vs infected with relaxed STAR filters

Same analysis as `../run_paired_samples/`, with STAR run with
`--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3`. The mapping test
(`../run_paired_samples/mapping_test/`) showed these settings recover reads lost to the high
divergence between our insects and the reference assembly.

- Pipeline: nf-core/rnaseq 3.21.0 (same version as the strict run), launched from `rnaseq/`
  with `-c relaxed_star.config` (see `rnaseq/run_nextflow.sh`); samplesheet
  `raw_reads/samples_paired_relaxed.csv`.
- Downstream (`star_salmon/deseq2_qc/`, conda env `R_popstat_jorge`): `bash run_downstream.sh`
  - `ruv_dea.r`: RUVs k=1 + DESeq2 (`~ W_1 + group`, lfcThreshold 1, greaterAbs, control vs infected)
  - `go_enrichment.r`: topGO BP classic Fisher, BH
  - `kegg_enrichment.r`: clusterProfiler enrichKEGG + GO–KEGG networks (from `kegg_claude.r`)
  - `ruv_exploration.r`: PCA plots of the RUVr / RUVg / RUVs corrections

## Two issues found in the strict analysis

1. **GO p-value filter.** topGO's `GenTable` returns p-values as text, and `filter(Classic < 0.05)`
   compared them as text, so terms printed in scientific notation were dropped. These were the
   most significant ones. Strict down lost 14 terms (e.g. *defense response to bacterium*,
   p = 3e-10; *innate immune response*) and strict up lost 4. `go_enrichment.r` uses the numeric
   p-values.
2. **RUV input.** The strict `deseq2.dds.RData` comes from an earlier run: its 30,491 rows are
   `gene-`/`rna-` IDs, while the quantification has 14,137 `DIATSA_LOCUS` genes. RUVs k=1 was
   estimated from it. Estimated from the run's own counts (nf-core `deseq2_qc.r` re-run on its
   `salmon.merged.gene_counts_length_scaled.tsv`), W_1 separates replicate 1 from replicates 2–3,
   which is the unwanted variation the correction was meant to remove. The strict DEGs then go
   from 147/82 to 266/147, and the originals are kept (147/147 up, 80/82 down).

Both issues are now fixed in `../run_paired_samples/star_salmon/deseq2_qc/`. The first version is
kept in its `original_analysis/` folder.

## Folders

| Folder | Content |
|---|---|
| `star_salmon/deseq2_qc/` | **relaxed** DEA, GO, KEGG, networks |
| `../run_paired_samples/star_salmon/deseq2_qc/` | corrected strict analysis (**baseline for the STAR comparison**) |
| `strict_reanalysis_original_ruv/` | strict counts with the original RUV input (reproduces the first 147/82 DEGs), fixed GO, current KEGG |
| `strict_vs_relaxed_summary.md/.tsv` | relaxed vs the corrected strict analysis |
| `strict_original_ruv_vs_relaxed_summary.md/.tsv` | relaxed vs `strict_reanalysis_original_ruv` |

KEGG enrichment downloads the current KEGG database. Some pathway names changed since the strict
analysis (e.g. *Cell adhesion molecules* is now *Cell adhesion molecule (CAM) interaction*).
