#!/bin/bash
##
## This script filters low quality reads
##
source ../variable.py
##
SP=$sp
PATH=$path
##
###############################
## Construct the file clean.sh file:
less $PATH/$SP/raw_seq_dir.txt | while read a b c; do echo "SOAPnuke filter -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -1 $b -2 $c -G -Q 2 -l 10 -q 0.2 -E 60 -5 0 -M 2 -o $PATH/$SP/trimmed_seq/$(basename $b '_read_1.fq.gz') -C $(basename $b '.fq.gz').clean.fq.gz -D $(basename $c '.fq.gz').clean.fq.gz"; done > clean.sh
##
## Split the files to submit them:
bash submit_split.sh
##
chmod u+wrx clean_splitted.sh
bash clean_splitted.sh
##
## Now the clean files are created in $PATH/£SP/trimmed_seq
