#!/bin/bash
##
## This script evaluate the insert size of the PE sequences by:
## sub sampling each reads
## Use perl 5.22.0
source ../variable.py
SP=$sp
RAW_DIR=raw_seq_dir_w.txt
##
## 0.smp sampling: ###################################################################################################################################################
## Import your directories from windows in 'raw_seq_dir_w.txt'
## Set up compatibility
tr -d '\15\32' < $PATH/$SP/$RAW_DIR > $PATH/$SP/raw_seq_dir.txt
## Create directories:
mkdir -p 0.smp/smp_seq/
## Construct the sampling.sh file:
echo "#!/bin/bash" > 0.smp/sampling.sh
echo "#SBATCH --partition normal" >> 0.smp/sampling.sh
echo "#SBATCH --mem-per-cpu 1G" >> 0.smp/sampling.sh
## Reads 1:
cut -f 2 $PATH/$SP/raw_seq_dir.txt | while read p; do echo "perl sampleReads.pl $p | gzip -c > 0.smp/smp_seq/$(basename $p '.fq.gz').smp.fq.gz"; done >> 0.smp/sampling.sh
## Reads 2:
cut -f 3 $PATH/$SP/raw_seq_dir.txt | while read p; do echo "perl sampleReads.pl $p | gzip -c > 0.smp/smp_seq/$(basename $p '.fq.gz').smp.fq.gz"; done >> 0.smp/sampling.sh
##
## Allow the execution of the file:
chmod u+x 0.smp/sampling.sh
##
## Submit the file:
sbatch 0.smp/sampling.sh
##
## This created a .smp.fa.gz file for all the lanes both ends with 1,000,000 reads
