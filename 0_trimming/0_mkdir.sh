## This script creates all the directories needed for the pipeline.
SPECIES=Macaca_mulatta
PATH=your_path

############################################################################################################
mkdir $PATH/$SPECIES/

mkdir $PATH/$SPECIES/ref_fasta

mkdir $PATH/$SPECIES/fastq_files
mkdir $PATH/$SPECIES/fastq_files/catenate

mkdir $PATH/$SPECIES/trimmed_seq

mkdir $PATH/$SPECIES/bam_files
mkdir $PATH/$SPECIES/bam_files/coverage
mkdir $PATH/$SPECIES/bam_files/coverage/cov.log
mkdir $PATH/$SPECIES/bam_files/inter_bam
mkdir $PATH/$SPECIES/bam_files/map.log
mkdir $PATH/$SPECIES/bam_files/merge.log
mkdir $PATH/$SPECIES/bam_files/recal.log
mkdir $PATH/$SPECIES/bam_files/summary
mkdir $PATH/$SPECIES/bam_files/uniq_rmdup.log

mkdir $PATH/$SPECIES/vcf_files
mkdir $PATH/$SPECIES/vcf_files/inter_vcf
mkdir $PATH/$SPECIES/vcf_files/back_com.log
mkdir $PATH/$SPECIES/vcf_files/call.log
mkdir $PATH/$SPECIES/vcf_files/combine.log
mkdir $PATH/$SPECIES/vcf_files/gather.log

mkdir $PATH/$SPECIES/vcf_handling
mkdir $PATH/$SPECIES/vcf_handling/log_file
mkdir $PATH/$SPECIES/vcf_handling/inter_vcf

mkdir $PATH/$SPECIES/de_novo_mutation
mkdir $PATH/$SPECIES/de_novo_mutation/log_file
mkdir $PATH/$SPECIES/de_novo_mutation/inter_vcf
