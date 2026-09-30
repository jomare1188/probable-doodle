

source ~/.bashrc
conda activate nextflow_25.10.4
#run 1
#nextflow run nf-core/rnaseq --input ../raw_reads/samples.csv --outdir run_paired_samples --gtf ../reference_genomes/diatraea_saccharalis/genomic.gtf --fasta ../reference_genomes/diatraea_saccharalis/GCA_918026875.4_PGI_DIATSA_v4_genomic.fna.gz --aligner star_salmon --skip_qc false -resume -profile conda


# run fusarium
nextflow run nf-core/rnaseq --input ../raw_reads/little_fungi_metadata.csv --outdir little_fungi --gtf ../reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.gtf.gz --fasta ../reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.fa.g
z --aligner star_salmon --skip_qc false -resume -profile conda

#run sugarcane
nextflow run nf-core/rnaseq --input /home/diegoj/rnaseq_diatraea/raw_reads/samples_fusarium_vs_sugarcane.csv --outdir run_fusarium_vs_sugarcane --gtf /home/diegoj/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf --fasta /home/diegoj/rnaseq_diatraea/reference_genomes/sugarcane/assembly/SofficinarumxspontaneumR570_771_v2.0.fa.gz --aligner star_salmon --skip_qc false -resume -profile conda -process.maxForks=3 --extra_star_align_args "--limitSjdbInsertNsj 10000139836" 

# next run will not need -process.maxForks, now star process has its own label and run only one star alignment at once utilizing all the server

# nextflow sugarcane too

nextflow run nf-core/rnaseq --input /home/diegoj/rnaseq_diatraea/raw_reads/samples_mock_cane_vs_fusarium-diatrea.csv --outdir run_mock-cane_vs_diatrea-fusarium --gtf /home/diegoj/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf --fasta /home/diegoj/rnaseq_diatraea/reference_genomes/sugarcane/assembly/SofficinarumxspontaneumR570_771_v2.0.fa.gz --aligner star_salmon --skip_qc false -resume -profile conda --extra_star_align_args "--limitSjdbInsertNsj 100139836"


# more sugarcane
samples_sugarcane-mock_diatrea.csv

nextflow run nf-core/rnaseq --input /home/diegoj/rnaseq_diatraea/raw_reads/samples_sugarcane-mock_diatrea.csv --outdir run_mock-cane_vs_diatrea --gtf /home/diegoj/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf --fasta /home/diegoj/rnaseq_diatraea/reference_genomes/sugarcane/assembly/SofficinarumxspontaneumR570_771_v2.0.fa.gz --aligner star_salmon --skip_qc false -resume -profile conda --extra_star_align_args "--limitSjdbInsertNsj 100139836"

# more sugarcane t3 vs t4 LEAF

nextflow run nf-core/rnaseq --input /dados04/jorge/rnaseq_diatraea/raw_reads/samples_leaf_t3_vs_t4.csv --outdir t3_vs_t4/ --gtf /dados04/jorge/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf --fasta /dados04/jorge/rnaseq_diatraea/reference_genomes/sugarcane/assembly/SofficinarumxspontaneumR570_771_v2.0.fa.gz --aligner star_salmon --skip_qc false -resume -profile conda --extra_star_align_args "--limitSjdbInsertNsj 100139836" -c star_align.config

# more sugarcane t3 vs t4 STEM
nextflow run nf-core/rnaseq --input /dados04/jorge/rnaseq_diatraea/raw_reads/samples_stem_t3_vs_t4.csv --outdir stem_t3_vs_t4/ --gtf /dados04/jorge/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf --fasta /dados04/jorge/rnaseq_diatraea/reference_genomes/sugarcane/assembly/SofficinarumxspontaneumR570_771_v2.0.fa.gz --aligner star_salmon --skip_qc false -resume -profile conda --extra_star_align_args "--limitSjdbInsertNsj 100139836" -c star_align.config -c conservative.config

# Diatraea paired samples again with relaxed STAR filters (mapping_test mode B: ~80% mapped instead of ~42%)
# nf-core/rnaseq 3.21.0 (same version as run_paired_samples) with moderate resources.
# Not the local clone in this folder: its edited modules/nf-core/star/align/main.nf (tabs in the
# versions heredoc) writes an invalid versions.yml and the run fails after STAR_ALIGN.
nextflow run nf-core/rnaseq -r 3.21.0 --input ../raw_reads/samples_paired_relaxed.csv --outdir run_paired_samples_relaxed_star --gtf ../reference_genomes/diatraea_saccharalis/genomic.gtf --fasta ../reference_genomes/diatraea_saccharalis/GCA_918026875.4_PGI_DIATSA_v4_genomic.fna.gz --aligner star_salmon --skip_qc false -resume -profile conda -c relaxed_star.config --extra_star_align_args '--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3'
