# Strict vs relaxed STAR: Diatraea control vs infected

Strict = default nf-core STAR filters; relaxed = `--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3`. Both analysed with the same scripts (RUVs k=1 DESeq2, topGO with numeric p-values, current KEGG).

## Mapping (per sample)

| Sample | Input pairs | STAR total mapped % strict | relaxed | Unique % strict | relaxed | Salmon reads strict (M) | relaxed (M) | Salmon gain |
|---|---|---|---|---|---|---|---|---|
| control_rep1 | 11.8 M | 44.2 | **80.0** | 42.1 | 72.6 | 3.96 | 5.55 | 1.40x |
| control_rep2 | 23.7 M | 42.0 | **78.9** | 40.0 | 71.0 | 7.43 | 10.55 | 1.42x |
| control_rep3 | 25.3 M | 40.9 | **76.0** | 38.9 | 68.7 | 7.71 | 10.98 | 1.42x |
| infected_rep1 | 12.6 M | 46.4 | **81.6** | 44.2 | 74.1 | 4.39 | 5.95 | 1.35x |
| infected_rep2 | 25.4 M | 44.1 | **80.6** | 42.1 | 72.8 | 8.32 | 11.37 | 1.37x |
| infected_rep3 | 27.3 M | 42.7 | **76.8** | 40.7 | 69.8 | 8.65 | 11.84 | 1.37x |

## Differentially expressed genes (RUVs k=1, |log2FC| > 1, padj < 0.05)

| Set | Strict | Relaxed | Shared | Strict only | Relaxed only |
|---|---|---|---|---|---|
| up (higher in control) | 147 | 331 | 139 | 8 | 192 |
| down (higher in infected) | 82 | 143 | 63 | 19 | 80 |

Genes tested: strict 10352, relaxed 11084, shared 10286. Pearson r of log2FC over shared genes: **0.885**.

## Enrichment

| Analysis | Strict | Relaxed | Shared | Strict only | Relaxed only |
|---|---|---|---|---|---|
| GO up | 200 | 243 | 119 | 81 | 124 |
| GO down | 64 | 87 | 36 | 28 | 51 |
| KEGG up | 5 | 2 | 2 | 3 | 0 |
| KEGG down | 5 | 4 | 2 | 3 | 2 |
