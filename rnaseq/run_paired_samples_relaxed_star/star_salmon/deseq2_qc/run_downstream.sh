#!/usr/bin/env bash
# Downstream analysis of an nf-core/rnaseq star_salmon run (Diatraea control vs infected):
#   1. RUVs (k = 1) + DESeq2 DEA   -> up_regulated.csv / down_regulated.csv
#   2. topGO GO enrichment         -> GO_up/down.csv, .svg, .pdf
#   3. KEGG enrichment + networks  -> kegg_up/down.*, gene_network_up/down.*
#   4. PCA plots of the RUVr / RUVg / RUVs corrections that were tried
# Run from <run>/star_salmon/deseq2_qc/ (defaults use paths relative to it).
#
# Usage: bash run_downstream.sh [samplesheet]
set -euo pipefail

SAMPLESHEET=${1:-../../../../raw_reads/samples_paired_relaxed.csv}
ROOT=../../../..
RSCRIPT=/home/genomics/miniconda3/envs/R_popstat_jorge/bin/Rscript
cd "$(dirname "$0")"

$RSCRIPT ruv_dea.r .. $SAMPLESHEET .
$RSCRIPT go_enrichment.r . $ROOT/panzzer/annot_01/gene_go_annotations.txt .
$RSCRIPT kegg_enrichment.r . . $ROOT
$RSCRIPT ruv_exploration.r .. $SAMPLESHEET .
