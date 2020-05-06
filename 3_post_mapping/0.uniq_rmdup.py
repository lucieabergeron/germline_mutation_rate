# -*- coding: utf-8 -*-
"""
This script keeps only the uniq mapped reads,
and remove the duplicates.
It finish with the indexing of this final file.

"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
directory = "{}/{}/bam_files/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}_sorted.merged.addg.bam".format(path, sp, name)


# The function:
def uniq_rmdup(input_bam, output_u, output_d, direct, sp):
    """The samtools keep uniq function."""
    uniq_cmd = "samtools view -h {} | ".format(input_bam)
    uniq_cmd += "grep -v -e 'XA:Z:' -e 'SA:Z:' | "
    uniq_cmd += "samtools view -b > {}".format(output_u)
    """The remove duplicate function"""
    rmdup_cmd = "gatk --java-options \"-XX:ParallelGCThreads=16 -Xmx100g -Djava.io.tmpdir=/scratch/$SLURM_JOBID/\" MarkDuplicates "
    rmdup_cmd += "-I= {} ".format(output_u)
    rmdup_cmd += "-OUTPUT= {} ".format(output_d)
    rmdup_cmd += "-REMOVE_DUPLICATES=true "
    rmdup_cmd += "-METRICS_FILE=metrics.txt "
    rmdup_cmd += "-MAX_FILE_HANDLES=800 "
    rmdup_cmd += "-TMP_DIR=/scratch/$SLURM_JOBID/ "
    rmdup_cmd += "-READ_NAME_REGEX=null"
    """The samtools index function."""
    index_cmd = "samtools index {}".format(output_d)
    """Create a .sh files with the keep uniq and remove duplicates functions."""
    file = open('{}_uniq_rmdup_idx.sh'.format(direct),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 128G \n')
    file.write('#SBATCH -c 16 \n')
    file.write('#SBATCH --time=30:00:00 \n')
    file.write(uniq_cmd)
    file.write('\n')
    file.write('echo \" Keeping unique IS DONE ##########################################################\" \n')
    file.write(rmdup_cmd)
    file.write('\n')
    file.write(index_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}_uniq_rmdup_idx.out {}_uniq_rmdup_idx.sh".format(direct, direct)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# For all merged bam files: keep uniq reads and remove duplicates.
for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}_sorted.merged.addg.uniq.bam".format(directory, name)):
        print("The uniq file for {} already exists".format(name))
        if os.path.exists("{}{}_sorted.merged.addg.uniq.rmdup.bam".format(directory, name)):
            print("The rmdup file for {} also exists ---> DONE".format(name))
        else:
            print("BUT the rmdup file for {} doesn't exists there is a problem".format(name))
    else:
        print("The uniq file for {} doesn't exist --> submit the function".format(name))
        uniq_rmdup(input_bam="{}{}_sorted.merged.addg.bam".format(directory, name), output_u="{}{}_sorted.merged.addg.uniq.bam".format(directory, name), output_d= "{}{}_sorted.merged.addg.uniq.rmdup.bam".format(directory, name), direct="{}{}".format(directory, name), sp=sp)
