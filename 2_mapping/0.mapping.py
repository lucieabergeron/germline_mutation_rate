# -*- coding: utf-8 -*-
"""
This script send the mapping jobs and save a txt files with the directories of the mapped sequences
"""
##################################################
# What you need ##################################
##################################################
# Packages:
import subprocess
import os
import pandas as pd
#import datetime
from variable import *

# Variables:
refDir="{}/{}/ref_fasta".format(path,sp)
print('\n The reference genome is located in: {}/{}.fa \n'.format(refDir, refGenome))

# Dictionary made of tuples of all read_1, read_2 for a given individual
f = open('{}/{}/clean_seq_dir.txt'.format(path, sp))
name2fastq_pairs = {}
for line in f:
    name, fq1, fq2 = line.split()
    if name.lower() == 'individual':
        continue
    if name not in name2fastq_pairs:
        name2fastq_pairs[name] = []
    name2fastq_pairs[name].append((fq1, fq2))

# The insert size of each lane:
ins_size = pd.read_csv('2.eva/res_ins_fin.txt', sep=' ', index_col=0)

# The function:
def submit_bwa_map(insert_size, ref_genome, read_1, read_2, output):
    """The `bwa mem` function."""
    map_cmd = "bwa mem -t 24 -I {} {} {} {}".format(insert_size, ref_genome, read_1, read_2)
    map_cmd += " | samtools sort -m 5G -@24 -O bam -T {} -o {}.bam".format(output, output)
    """Create a .sh files with the `bwa mem` function."""
    file = open('{}.sh'.format(output),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 256G \n')
    file.write('#SBATCH -c 24 \n')
    file.write('#SBATCH --time=12:00:00 \n')
    file.write(map_cmd)
    file.write('\n')
    file.close()
    """Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}.out {}.sh".format(output, output)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################
# For each sample, you map all the lanes and give them sample_name_nb_of_lane
# A file with all the output mapped lane for merging.py
all_lane=[]
bam_files_directories = open("{}/{}/bam_files_directories.txt".format(path,sp), "w")
for name in name2fastq_pairs:
    bam_dir = [name]
    for i in range(len(name2fastq_pairs[name])):
        read_1, read_2 = name2fastq_pairs[name][i]
        seq = os.path.basename(read_1)[:-19]
        ##
        ins_mean=round(ins_size.loc[seq]['MEAN'], 2)
        ins_sd=round(ins_size.loc[seq]['SD'], 2)
        ins_max=int(round(ins_size.loc[seq]['MAX']))
        ins_min=int(round(ins_size.loc[seq]['MIN']))
        bamfile_dir = "{}/{}/bam_files/{}_{}.sorted".format(path, sp, name, i)
        if os.path.exists(bamfile_dir + '.bam'):
            print("Mapping not sent because " + bamfile_dir + ".bam exist already")
        else:
            print(seq + " --> " + bamfile_dir + '.bam')
            all_lane.append(seq)
            submit_bwa_map(insert_size='{},{},{},{}'.format(ins_mean, ins_sd, ins_max, ins_min),ref_genome='{}/{}'.format(refDir, refGenome), read_1=read_1, read_2=read_2, output=bamfile_dir)
        bam_dir.append(bamfile_dir + '.bam')
    bam_dir = ' '.join(bam_dir)
    bam_files_directories.write(bam_dir + "\n")

all_lane = ' '.join(all_lane)
print("This lanes have been sent for mapping: {}".format(all_lane))
bam_files_directories.close()
print("The directories of the mapped files are stored in {}/{}/bam_files_directories.txt".format(path, sp))
