# -*- coding: utf-8 -*-
"""
This script catenate the GVCF files per chromosome for all individuals.
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *
import pandas as pd

# Directories:
direct = "{}/{}/vcf_files/".format(path, sp)
chrom_dir = "{}/{}/".format(path, sp)


# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}_sorted.merged.addg.uniq.rmdup.bam".format(path, sp, name)

# And VCF
f = open('{}/{}/vcf_files.txt'.format(path, sp))
vcf_dir = {}
for line in f:
    name = line.split()[0]
    if name not in vcf_dir:
        vcf_dir[name] = []
    vcf_dir[name] = "{}/{}/vcf_files/{}".format(path, sp, name)


# Import chrom names:
chrom_name = pd.read_csv('{}chromosomes.txt'.format(chrom_dir),sep=' ', index_col=None, header=None)

# The function:
def combine(chrom_old, chrom_new, direct):
    """Combine all samples with GenomicsDBImport"""
    combine_cmd = "gatk --java-options \"-XX:ParallelGCThreads=1 -Xmx100g -Djava.io.tmpdir=/scratch/$SLURM_JOBID/\" GenomicsDBImport "
    for i in vcf_dir:
        combine_cmd += "--variant {}{} ".format(direct,i)
    combine_cmd += "--TMP_DIR /scratch/$SLURM_JOBID/ "
    combine_cmd += "--genomicsdb-workspace-path /scratch/$SLURM_JOBID/genomicDBI_{} ".format(chrom_new)
    combine_cmd += "-L {} ".format(chrom_old)
    """Create a .sh files with the combine variant functions."""
    file = open('{}combine_genomicDBImport_{}.sh'.format(direct, chrom_new),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition short,normal \n')
    file.write('#SBATCH --mem 110G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=10:00:00 \n')
##    file.write('#SBATCH --time=250:00:00 \n')
    file.write(combine_cmd)
    file.write('\n')
    file.write('cp -a /scratch/$SLURM_JOBID/genomicDBI_{} {}genomicDBI_{}'.format(chrom_new, direct, chrom_new))
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}combine_genomicDBImport_{}.out {}combine_genomicDBImport_{}.sh".format(direct, chrom_new, direct, chrom_new)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# For each chromosome one function:
list_exist=[]
for file in list(vcf_dir.values()):
    list_exist.append(os.path.exists(file))
if all(list_exist):
    print("\t All the res.g.vcf files exist --> call variant done")
    mv_call = "mv {}*_call_res_g* {}call.log".format(direct, direct)
    subprocess.call(mv_call, shell=True)
    print("\t Move the call log files")
    for line in range(0, nb_chrom):
        chrom_old=chrom_name.loc[line,0]
        chrom_new=chrom_name.loc[line,1]
        vcf_dir=[]
        for name in bamfile_dir:
            vcf_dir.append("{}_{}_res.g.vcf".format(name, chrom_new))
        print("Combine for chromosome {} called {}".format(chrom_new, chrom_old))
        combine(chrom_old=chrom_old, chrom_new=chrom_new, direct=direct)
else:
    print("\t Some res.g.vcf file are missing --> PROBLEM")


