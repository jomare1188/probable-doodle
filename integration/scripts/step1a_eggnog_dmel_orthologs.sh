#!/usr/bin/env bash
# Step 1a - Drosophila melanogaster orthologs for the Diatraea saccharalis proteome.
# The published annotation (eggnog/annotation/proteins.emapper.emapper.annotations; emapper 2.1.3,
# eggNOG DB 5.0.2, DIAMOND mode, default settings, proteome GCA_918026875.4 PGI_DIATSA_v4) is kept
# as is for KOs and gene names. Here its seed-ortholog hits are re-annotated WITHOUT a new search
# (-m no_search) to list the orthologs of each query in D. melanogaster (NCBI taxid 7227).
# Run from the repository root with the eggnog_env env:  bash scripts/step1a_eggnog_dmel_orthologs.sh
# Requires emapper.py on PATH and EGGNOG_DATA_DIR (config.sh): the eggNOG 5.0.2 data folder made by
# download_eggnog_data.py (eggnog.db, eggnog_proteins.dmnd, and eggnog.taxa.db or eggnog.taxa.tar).
set -euo pipefail

: "${EGGNOG_DATA_DIR:?set EGGNOG_DATA_DIR in config.sh}"
SEEDS=data/annotation/proteins.emapper.emapper.seed_orthologs
OUT=results/step1_gene_annotation
mkdir -p "$OUT"

# eggNOG data dir for this step: eggnog.db and the DIAMOND db are symlinked, and the NCBI
# taxonomy db (needed for --target_taxa) is linked or unpacked from eggnog.taxa.tar here, so the
# shared eggNOG folder is not modified.
DATA=$OUT/eggnog_data
mkdir -p "$DATA"
for f in eggnog.db eggnog_proteins.dmnd; do ln -sfn "$EGGNOG_DATA_DIR/$f" "$DATA/$f"; done
if [ ! -f "$DATA/eggnog.taxa.db" ]; then
  if [ -f "$EGGNOG_DATA_DIR/eggnog.taxa.db" ]; then
    ln -sfn "$EGGNOG_DATA_DIR/eggnog.taxa.db" "$DATA/eggnog.taxa.db"
    ln -sfn "$EGGNOG_DATA_DIR/eggnog.taxa.db.traverse.pkl" "$DATA/eggnog.taxa.db.traverse.pkl"
  else
    tar -xf "$EGGNOG_DATA_DIR/eggnog.taxa.tar" -C "$DATA"
  fi
fi

emapper.py -m no_search \
  --annotate_hits_table "$SEEDS" \
  --data_dir "$DATA" \
  --report_orthologs --target_orthologs all --target_taxa 7227 \
  --output dmel_orthologs --output_dir "$OUT" --cpu "${THREADS:-16}" --override
# compressed for storage; step 1b reads the .gz file
gzip -f "$OUT"/dmel_orthologs.emapper.orthologs "$OUT"/dmel_orthologs.emapper.annotations
