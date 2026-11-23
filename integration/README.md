# *Diatraea saccharalis* multi-omics integration

Knowledge-based, summary-level integration of three omics measured in *Diatraea saccharalis* larvae fed sugarcane with (infected) or without (control) *Fusarium verticillioides*:

- **transcriptome**: published DESeq2 results;
- **metabolome**: 18 published metabolites;
- **saliva microbiome**: Oxford Nanopore amplicon reads, re-classified here.

Each omic comes from different larvae, so no sample-level method is used (no correlation, DIABLO, MOFA…). The omics are linked only through prior knowledge: KEGG reactions and pathways, MicrobiomeAnalyst, MetOrigin, and the literature. Every effect is expressed as **infected vs control**.

## Folder layout

| folder | content |
|---|---|
| `scripts/` | all analysis scripts, numbered in run order |
| `data/` | input data (published results and annotation), described below |
| `resources/` | frozen copies of web resources: PubChem / KEGG REST responses, MicrobiomeAnalyst reference libraries, SILVA taxonomy build script |
| `web_inputs/` | lists uploaded to the web tools |
| `web_results/` | exports from the MicrobiomeAnalyst and MetOrigin websites |
| `results/` | outputs of every step, one folder per step, each with `sessionInfo.txt` for the R steps |
| `deliverables/` | Excel workbook with all tables (README sheet + 25 sheets), final figures, Cytoscape node / edge tables |
| `envs/` | conda environments and Python requirements |

## Pipeline

| order | script | step | main outputs |
|---|---|---|---|
| 1 | `step0_qc_harmonize.R` | 0 QC and harmonisation; published DESeq2 log2FC (control/infected) flipped to infected vs control | `results/step0/` |
| 2 | `step0a_microbiome_cascade.py` | 0a hierarchical kraken2 classification, each read counted once (SILVA 138.2 → NCBI 16S/18S/28S/ITS → PlusPFP-16) | `results/step0a_microbiome_cascade/` genus tables, per-read assignment |
| 3 | `step0a_compare_tables.py` | 0a comparison with the earlier published genus tables (descriptive) | `comparison_with_published_tables.tsv` |
| 4 | `step0b_microbiome_DA.R` | 0b differential abundance: ALDEx2 (primary), ANCOM-BC2, strict cascade variant, PERMANOVA (descriptive) | `results/step0b_microbiome_DA/` |
| 5 | `step1a_eggnog_dmel_orthologs.sh` | 1a *D. melanogaster* orthologs from the published eggNOG seed hits (`-m no_search`) | `dmel_orthologs.emapper.orthologs.gz` |
| 6 | `step1b_gene_annotation.R` | 1b gene-level KO / name / fly-ortholog table and coverage | `results/step1_gene_annotation/` |
| 7 | `step2_metabolite_ids.py`, `step2b_finalize_metabolite_ids.py` | 2 metabolite IDs (KEGG, HMDB, ChEBI, PubChem), two independent routes plus curated overrides | `results/step2_metabolite_ids/` |
| 8 | `step3_multigsea.R` | 3 multiGSEA on KEGG fly pathways (genes + metabolites), Stouffer combination | `results/step3_multigsea/` |
| 9 | `step0c_web_inputs.py` | 0c taxon and metabolite lists for the web tools | `web_inputs/` |
| — | MicrobiomeAnalyst, MetOrigin (websites) | manual runs with the files in `web_inputs/` | `web_results/` |
| 10 | `step4a_microbiomeanalyst_local.R` | 4a MicrobiomeAnalyst list-input analyses reproduced locally from the server's reference libraries | `results/step4a_microbiomeanalyst_local/` |
| 11 | `step4b_kegg_reaction_tables.py` | 4b KEGG compound → reaction → KO and compound → pathway tables | `results/step4b_kegg/` |
| 12 | `step4c_integration.R` | 4 integration: (a) metabolite evidence, (b) pathway × omics matrix, (c) direction agreement, (d) knowledge network, (e) literature links | `results/step4_integration/` |
| 13 | `step5a_validate_outputs.R` | independent recomputation of key results of steps 0–3 (stops on any FAIL) | `results/validation_steps0-3.csv` |
| 14 | `step5b_deliverables.R` | workbook, figures, Cytoscape tables | `deliverables/` |

## Input data

| file | origin |
|---|---|
| `data/transcriptome/up_regulated.csv`, `down_regulated.csv` | published DESeq2 DE lists (padj < 0.05), nf-core/rnaseq STAR–Salmon run with relaxed STAR settings, model `~ W_1 + group` (RUVs) |
| `data/transcriptome/dea_all_genes.csv` | the full DESeq2 table of the same run (all tested genes) |
| `data/metabolome/metabolome.xlsx` | published metabolite abundances, 6 control + 6 infected |
| `data/metabolome/statistical_results.csv` | published per-metabolite statistics (p, p_adj, log2FC) |
| `data/microbiome/abundance_table_ready.tsv`, `abundance_table_genus_CONSOLIDATED.tsv` | earlier genus tables, used only for the comparison in step 0a (the cascade table replaces them) |
| `data/annotation/proteins.emapper.emapper.annotations.gz`, `.seed_orthologs` | published eggNOG-mapper 2.1.3 annotation of the *D. saccharalis* proteome (GCA_918026875.4, PGI_DIATSA_v4; eggNOG 5.0.2, DIAMOND mode) |
| `data/annotation/gene2protein.txt` | gene ↔ protein table of the same genome annotation |

The raw saliva reads (6 samples) are not stored here because of their size (485 MB). Step 0a reads them from `RAW_READS_DIR`, which holds one sub-folder per sample (`control_1` … `infected_3`) with the basecalled `fastq.gz` files.

## Frozen web resources (`resources/`)

| folder | content | retrieved |
|---|---|---|
| `web_cache/step2_pubchem_kegg/` | PubChem PUG REST and KEGG REST responses of step 2 | 2026-09-30 |
| `web_cache/step4b_kegg/` | KEGG REST responses of step 4b (KEGG release 2026-09-30) | 2026-10-01 |
| `microbiomeanalyst_libs/` | reference libraries of the MicrobiomeAnalyst server (`https://www.microbiomeanalyst.ca/resources/libs/mmp/`), MD5 sums in `MD5SUMS.txt` | 2026-10-01 |
| `kraken2_silva_138.2/build_silva_taxonomy.pl` | Kraken2 script used to turn `tax_slv_ssu_138.2.txt` into the SILVA `names.dmp` / `nodes.dmp` | — |

With these caches in place no step queries PubChem or KEGG, so the outputs do not depend on later database updates. If a cache file is deleted, the script downloads it again. Step 3 uses the KEGG pathway data that graphite 1.52.0 distributes; the retrieval date of that data (2024-10-22) is printed in `results/step3_multigsea/sessionInfo.txt`.

## Software versions

R 4.4.3, ALDEx2 1.38.0, ANCOMBC 2.8.0, phyloseq 1.50.0, vegan 2.7.5, multiGSEA 1.16.0, fgsea 1.32.2, graphite 1.52.0, org.Dm.eg.db 3.20.0, qs 0.27.3, ComplexHeatmap 2.22.0, UpSetR 1.4.1, igraph 2.3.3, ggraph 2.2.2, openxlsx 4.2.9; kraken2 2.17.2, SeqKit 2.14.0, TaxonKit 0.20.0; eggNOG-mapper 2.1.3 (eggNOG 5.0.2, DIAMOND 2.0.15); Python 3.13.11, pandas 3.0.2. The complete package lists are in `envs/`.
