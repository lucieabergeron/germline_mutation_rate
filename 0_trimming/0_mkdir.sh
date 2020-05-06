## This script creates all the directories needed for the pipeline.
SPECIES=Macaca_mulatta

############################################################################################################
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/

mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/ref_fasta

mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/fastq_files
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/fastq_files/catenate

mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/trimmed_seq

mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/coverage
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/coverage/cov.log
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/inter_bam
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/map.log
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/merge.log
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/recal.log
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/summary
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/bam_files/uniq_rmdup.log

mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_files
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_files/inter_vcf
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_files/back_com.log
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_files/call.log
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_files/combine.log
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_files/gather.log

mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_handling
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_handling/log_file
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/vcf_handling/inter_vcf

mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/de_novo_mutation
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/de_novo_mutation/log_file
mkdir /home/lucie/MammalianMutation/faststorage/$SPECIES/de_novo_mutation/inter_vcf
