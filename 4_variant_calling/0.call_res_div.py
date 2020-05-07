# -*- coding: utf-8 -*-
"""
This script call variants in BR RESOLUTION for all individuals per chromosomes
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
ref_dir = "{}/{}/ref_fasta/{}.fa".format(path, sp, refGenome)
bam_dir = "{}/{}/bam_files/".format(path, sp)
vcf_dir = "{}/{}/vcf_files/".format(path, sp)
chrom_dir = "{}/{}/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}_sorted.merged.addg.uniq.rmdup.bam".format(path, sp, name)


# Import chrom names:
chrom_name = pd.read_csv('{}chromosomes.txt'.format(chrom_dir),sep=' ', index_col=None, header=None)


# The function:
def call_var(ref, in_bam, out_vcf, chrom_old, out_dir, chrom_new):
    """Happlotype caller function to call variants for each samples"""
    call_cmd = "gatk --java-options \"-XX:ParallelGCThreads=1 -Xmx90g -Djava.io.tmpdir=/scratch/$SLURM_JOBID/\" HaplotypeCaller "
    call_cmd += "-R= {} ".format(ref)
    call_cmd += "-I= {} ".format(in_bam)
    call_cmd += "-O= {} ".format(out_vcf)
    call_cmd += "-ERC BP_RESOLUTION "
    call_cmd += "-L {} ".format(chrom_old)
    call_cmd += "--dont-use-soft-clipped-bases "
    call_cmd += "--native-pair-hmm-threads 1 "
    call_cmd += "--TMP_DIR=/scratch/$SLURM_JOBID/ "
    """Create a .sh files with the calling variant functions."""
    file = open('{}_call_res_g_{}.sh'.format(out_dir, chrom_new),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 100G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=40:00:00 \n')
    file.write(call_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}_call_res_g_{}.out {}_call_res_g_{}.sh".format(out_dir, chrom_new, out_dir, chrom_new)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# For all merged bam files: keep uniq reads and remove duplicates.

vcf_files_dir = open("{}/{}/vcf_files.txt".format(path, sp), "w")
for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}_res.g.vcf".format(vcf_dir, name)):
        print("\t The res.g.vcf file for {} already exists --> CALL VARIANT DONE".format(name))
    else:
        print("\t The res.g.vcf file for {} doesn't exist --> submit the function".format(name))
        for line in range(0, nb_chrom):
            chrom_old=chrom_name.loc[line,0]
            chrom_new=chrom_name.loc[line,1]
            call_var(ref=ref_dir, in_bam=bamfile_dir[name], out_vcf="{}{}_{}_res.g.vcf".format(vcf_dir, name, chrom_new), chrom_old=chrom_old, out_dir="{}{}".format(vcf_dir, name), chrom_new=chrom_new)
            if "/" in chrom_old:
                list_nb=chrom_old.split('/')[8]
                list_nb_1=list_nb.split('.')[0]
            else:
                list_nb_1=chrom_new
            vcf_files_dir.write(name + "_" + list_nb_1 +"_res.g.vcf \n")
vcf_files_dir.close()

