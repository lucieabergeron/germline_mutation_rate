#!/bin/bash
##
## This script runs the summary_SOP for all sequences
##
source ../variable.py
##
SP=$sp
PATH=$path
##
## Remove the output of cleanning:
rm submit_clean*.out
## Again build a bash with all the bash for that
ls -d $PATH/$SP/trimmed_seq/*/ | while read p; do echo "bash summary_SOAPnuke.sh $PATH/$SP/trimmed_seq/$(basename "$p" '/') $PATH/$SP/trimmed_seq"; done > run_summary_SOAPnuke.sh
##
##
chmod u+rwx run_summary_SOAPnuke.sh
bash run_summary_SOAPnuke.sh
##
## Now you should have summary_SOAPnuke.txt in the output directory : $PATH/$SP/trimmed_seq/
##
## Write the directories of the trimmed files
less $PATH/$SP/raw_seq_dir.txt | while read a b c; do echo -e "$a\t$(echo $b | sed -r "s:.fq.gz:.clean.fq.gz:g; s:fastq_files:trimmed_seq/$(basename "$b" '_read_1.fq.gz'):g;")\t$(echo $c | sed -r "s:.fq.gz:.clean.fq.gz:g; s:fastq_files:trimmed_seq/$(basename "$b" '_read_1.fq.gz'):g;")";done > $PATH/$SP/clean_seq_dir.txt
