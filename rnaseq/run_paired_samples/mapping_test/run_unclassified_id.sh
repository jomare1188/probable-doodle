#!/usr/bin/env bash
# Identify reads left unmapped (relaxed STAR) and unclassified (Kraken2 PlusPFP)
#  1. sensitive blastn of 10k reads vs D. saccharalis genome / transcripts
#  2. DIAMOND blastx of reads vs UniProt Swiss-Prot
#  3. Trinity assembly of all unclassified pairs, salmon quantification,
#     contig annotation (blastn genome, DIAMOND Swiss-Prot, TransDecoder + Pfam)
#  4. read-weighted category summary
set -euo pipefail

W=/tmp/claude-1004/-dados04-jorge-rnaseq-diatraea/b1f93fcc-685d-438d-8e75-221eba76ed5d/scratchpad/maptest
REPO=/dados04/jorge/rnaseq_diatraea
REF=$REPO/reference_genomes/diatraea_saccharalis
SPROT=/dados04/jorge/databases/uniprot/uniprot_sprot.dmnd
PFAM=/dados04/jorge/databases/pfam/Pfam-A.hmm
ENVS=/home/genomics/miniconda3/envs
DIAMOND=$ENVS/diamond_2.2.3/bin/diamond
THREADS=32
SAMPLES="control_rep1 control_rep2 control_rep3 infected_rep1 infected_rep2 infected_rep3"
cd $W; mkdir -p uid uncl_fix; cd uid

# 0. Kraken2 writes mates as "0:N"/"1:N"; Trinity needs /1 /2
for f in ../uncl/*_[12].fq; do
  b=$(basename $f .fq); m=${b##*_}
  [ -s ../uncl_fix/$b.fq ] || awk -v m=$m 'NR%4==1{split($1,a," "); print a[1]"/"m; next}{print}' $f > ../uncl_fix/$b.fq
done

# 1. read-level blastn vs Diatraea genome and transcripts
if [ ! -s uncl10k.fa ]; then
  for s in $SAMPLES; do seqkit fq2fa ../uncl/${s}_1.fq | seqkit replace -p '^(\S+).*' -r "${s}|\$1"; done \
    | seqkit seq -m 80 | seqkit shuffle -s 42 | seqkit head -n 10000 > uncl10k.fa
fi
cp $REPO/rnaseq/run_paired_samples/mapping_test/classified_100_reads.fasta ctrl100.fa
cp $REPO/rnaseq/run_paired_samples/mapping_test/unclassified_100_reads.fasta user100.fa
[ -s dsac_rna.fa ] || zcat $REF/GCA_918026875.4_PGI_DIATSA_v4_rna_from_genomic.fna.gz > dsac_rna.fa
[ -s db_genome.nsq ] || [ -s db_genome.00.nsq ] || makeblastdb -in ../genome.fa -dbtype nucl -out db_genome
[ -s db_rna.nsq ] || makeblastdb -in dsac_rna.fa -dbtype nucl -out db_rna
F="6 qseqid sseqid pident length qlen evalue bitscore"
for q in uncl10k user100 ctrl100; do
  for db in genome rna; do
    [ -s ${q}_vs_${db}.tsv ] || blastn -task blastn -query $q.fa -db db_$db -evalue 1e-5 \
           -max_target_seqs 1 -max_hsps 1 -outfmt "$F" -num_threads $THREADS > ${q}_vs_${db}.tsv
  done
  # 2. read-level DIAMOND blastx vs Swiss-Prot
  [ -s ${q}_vs_sprot.tsv ] || $DIAMOND blastx --more-sensitive -e 1e-5 -k 1 -p $THREADS -d $SPROT \
           -q $q.fa -o ${q}_vs_sprot.tsv -f 6 qseqid pident length evalue bitscore stitle --quiet
done

# 3. Trinity assembly of all unclassified pairs
L=$(ls ../uncl_fix/*_1.fq | paste -sd,); R=$(ls ../uncl_fix/*_2.fq | paste -sd,)
[ -s trinity_uncl.Trinity.fasta ] || [ -s trinity_uncl/Trinity.fasta ] || \
  ( export PATH=$ENVS/trinity/bin:$PATH; Trinity --seqType fq --left $L --right $R --CPU $THREADS \
       --max_memory 100G --output trinity_uncl > trinity.log 2>&1 )
T=$( [ -s trinity_uncl/Trinity.fasta ] && echo trinity_uncl/Trinity.fasta || echo trinity_uncl.Trinity.fasta )
cp $T contigs.fa
$ENVS/trinity/bin/TrinityStats.pl contigs.fa > trinity_stats.txt

# 3a. quantify all unclassified pairs against the assembly
$ENVS/trinity/bin/salmon index -t contigs.fa -i salmon_idx -p $THREADS > /dev/null 2>&1
$ENVS/trinity/bin/salmon quant -i salmon_idx -l A -1 ${L//,/ } -2 ${R//,/ } -p $THREADS \
       --validateMappings -o salmon_quant > salmon.log 2>&1

# 3b. contig annotation
blastn -task blastn -query contigs.fa -db db_genome -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 \
       -outfmt "$F" -num_threads $THREADS > contigs_vs_genome.tsv
$DIAMOND blastx --more-sensitive -e 1e-5 -k 1 -p $THREADS -d $SPROT -q contigs.fa \
       -o contigs_vs_sprot.tsv -f 6 qseqid pident length evalue bitscore stitle --quiet
( export PATH=$ENVS/transdecoder/bin:$PATH; $ENVS/transdecoder/opt/transdecoder/util/TransDecoder.LongOrfs -t contigs.fa -m 100 > transdecoder.log 2>&1 )
( export PATH=$ENVS/transdecoder/bin:$PATH; $ENVS/transdecoder/opt/transdecoder/util/TransDecoder.Predict -t contigs.fa --single_best_only >> transdecoder.log 2>&1 )
$ENVS/transdecoder/bin/hmmsearch --cpu $THREADS -E 1e-5 --domtblout contigs_pfam.domtbl \
       $PFAM contigs.fa.transdecoder.pep > /dev/null

# 4. categories and outputs
python3 $W/categorize_contigs.py
echo DONE
