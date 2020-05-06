#!/bin/bash
##
## This script evaluate the insert size of the PE sequences by:
## 2. single end mapping the reads,
##
source ../variable.py
##################################################
SP=$sp
PATH=$path
SP_REF=$refGenome
## create directories:
mkdir -p 1.map/map_seq
## 1.map mapping: ####################################################################################################################################################
## Reconstruction of mapping.sh file by running "map.sh path/ref_genome" (no .fa):
##
bash map.sh $PATH/$SP/ref_fasta/$SP_REF
##
## Then you run a bash file to split the single map and submit them in the same time:
bash submit_split.sh
##
## Running the mapping_splitted.sh to submit all the job:
bash 1.map/mapping_splitted.sh
##
## This creates 2 bam file: read 1 and read 2
