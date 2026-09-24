#!/usr/bin/env bash
# Relaxed STAR mapping test for Diatraea samples (reviewer comment on 42% mapping rate)
set -euo pipefail

W=/tmp/claude-1004/-dados04-jorge-rnaseq-diatraea/b1f93fcc-685d-438d-8e75-221eba76ed5d/scratchpad/maptest
REPO=/dados04/jorge/rnaseq_diatraea
REF=$REPO/reference_genomes/diatraea_saccharalis
THREADS=32
NPAIRS=1000000
cd $W

IMG=$W/star_2.7.11b.sif
[ -f $IMG ] || singularity pull $IMG docker://quay.io/biocontainers/star:2.7.11b--h5ca1c30_5
STAR="singularity exec -B $W,$REPO $IMG STAR"

# 1. genome index (same GTF as the pipeline)
if [ ! -f index/SAindex ]; then
  mkdir -p index
  zcat $REF/GCA_918026875.4_PGI_DIATSA_v4_genomic.fna.gz > genome.fa
  $STAR --runMode genomeGenerate --runThreadN $THREADS --genomeDir index \
        --genomeFastaFiles genome.fa --sjdbGTFfile $REF/genomic.gtf --sjdbOverhang 99
fi

# 2. subsample first NPAIRS pairs per sample
declare -A MAP=( [control_rep1]=interaction1_rep1 [control_rep2]=interaction1_rep2 [control_rep3]=interaction1_rep3
                 [infected_rep1]=interaction2_rep1 [infected_rep2]=interaction2_rep2 [infected_rep3]=interaction2_rep3 )
mkdir -p reads out
for s in "${!MAP[@]}"; do
  for r in 1 2; do
    f=reads/${s}_R$r.fq
    [ -s $f ] || { zcat $REPO/raw_reads/${MAP[$s]}_R${r}_paired.fq.gz || true; } | head -n $((NPAIRS*4)) > $f
  done
done

# nf-core rnaseq default STAR args (without BAM output)
BASE="--genomeDir index --runThreadN $THREADS --outSAMtype None --twopassMode Basic \
      --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --runRNGseed 0"

# 3. three mapping modes
for s in "${!MAP[@]}"; do
  $STAR $BASE --readFilesIn reads/${s}_R1.fq reads/${s}_R2.fq --outFileNamePrefix out/${s}.A_default.
  $STAR $BASE --readFilesIn reads/${s}_R1.fq reads/${s}_R2.fq --outFileNamePrefix out/${s}.B_relaxed. \
        --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
  $STAR $BASE --readFilesIn reads/${s}_R1.fq --outFileNamePrefix out/${s}.C_R1only.
done

# 4. summary
{
  printf "sample\tmode\tunique%%\tmulti%%\ttotal_mapped%%\ttoo_short%%\tother%%\tmismatch%%\n"
  for f in out/*.Log.final.out; do
    b=$(basename $f .Log.final.out); s=${b%.*}; m=${b##*.}
    g() { grep "$1" $f | cut -d'|' -f2 | tr -d ' \t%'; }
    u=$(g "Uniquely mapped reads %"); mu=$(g "% of reads mapped to multiple loci"); tm=$(g "% of reads mapped to too many loci")
    printf "%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\n" $s $m $u $mu $(echo "$u+$mu+$tm" | bc) \
      $(g "% of reads unmapped: too short") $(g "% of reads unmapped: other") $(g "Mismatch rate per base")
  done
} > summary.tsv
column -t summary.tsv
echo DONE
