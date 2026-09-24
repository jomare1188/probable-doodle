#!/usr/bin/env bash
# Classify reads left unmapped by relaxed STAR mapping (mode B) with kraken2 PlusPFP
set -euo pipefail

W=/tmp/claude-1004/-dados04-jorge-rnaseq-diatraea/b1f93fcc-685d-438d-8e75-221eba76ed5d/scratchpad/maptest
REPO=/dados04/jorge/rnaseq_diatraea
DB=/dados04/jorge/databases/kraken2/k2_pluspfp_16_GB_20260626
THREADS=32
cd $W
STAR="singularity exec -B $W,$REPO star_2.7.11b.sif STAR"
mkdir -p unmapped kraken

for s in control_rep1 control_rep2 control_rep3 infected_rep1 infected_rep2 infected_rep3; do
  # 1. relaxed mapping, keep unmapped pairs
  $STAR --genomeDir index --runThreadN $THREADS --outSAMtype None --twopassMode Basic \
        --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --runRNGseed 0 \
        --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
        --outReadsUnmapped Fastx \
        --readFilesIn reads/${s}_R1.fq reads/${s}_R2.fq --outFileNamePrefix unmapped/${s}.
  # 2. kraken2 on unmapped pairs
  conda run -n kraken2 kraken2 --db $DB --threads $THREADS --paired \
        --report kraken/${s}.report --output /dev/null \
        unmapped/${s}.Unmapped.out.mate1 unmapped/${s}.Unmapped.out.mate2
done
echo DONE
