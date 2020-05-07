# Estimation of germline mutation rate

This pipeline estimate germline mutation rate from NGS data in pedigree samples (mother, father and offspring). This pipeline was built to run on any species but the example given here is for 19 trio of rhesus macaques (*Macaca mulatta*). 
The raw sequences can be found on NCBI under the project number PRJNA588178.

## Overview

The pipeline was design for a slurm job scheduling system on Linux. The workflow is composed of python script writing and submitting bash jobs to the cluster. Using another job scheduling system would require some changes in the scripts. 

A first step (0_settings) informe on the directories to create and file to write before running the pipeline.
Then the pipeline consists of 6 steps that should be run in order. In each step a README.txt explains how to run the different scripts.

The steps are:
1. Trimming the raw fatq files and estimating insert size.
2. Mapping each lanes to the reference genome (bwa mem) and merging all the lanes for a single individual.
3. Processing the bam file to remove PCR duplicate reads or reads mapping to multiple location (+ quality control on coverage).
4. Calling variant with GATK4 by calling per individuals with HaplotypeCaller in BP-RESOLUTION, combining all gVCF files with CombineGVCFs, joint genotyping with GenotypeGVCF.
5. Filtering the sites and detecting mendelian violation on the standard vcf file.  
6. Filtering *de novo* mutation candidates, estimating the number of callable sites on the bp-resolution vcf file and calculating mutation rate.

## Requirements 

Default software used to run the pipeline are:
- bwa 0.7.15
- samtools 0.1.18
- soapnuke 1.5.6
- python 3.7.3
- java 1.8.0_222
- gatk 4.0.7.0
- bcftools 1.9
- R 3.5.1
- bcftools 1.9

For some script the following version were used (specified in the step/README.txt):
- samtools 1.2
- samtools 1.4.1
- samtools 0.1.19
- picard.jar 2.7.1
