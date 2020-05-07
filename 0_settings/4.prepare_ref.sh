# The reference genome should be indexed with:
bwa index -p {ref_genome} -a bwtsw {ref_genome}.fa
java -jar picard.jar CreateSequenceDictionary R={ref_genome}.fa O={ref_genome}.dict
samtools faidx {ref_genome}.fa
