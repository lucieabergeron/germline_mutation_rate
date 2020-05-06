# -*- coding: utf-8 -*-
"""
This script check if the mapping has created the files,
set the summary of each mapped lane,
merged the lane from the same library,
and add a group name to the merged mapped bamfile.

"""
##################################################
# What you need ##################################
##################################################
# Packages:
import subprocess
import os
from variable import *

# Directory:
directory = "{}/{}/bam_files/".format(path, sp)

# Dictionary made of tuples of all mapped bamfile for a sample
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name, *bamfile = line.split()
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name].append(bamfile)


# The functions:
def merge_add(name, inputs, ID, output):
    """The samtools merge function."""
    mg_cmd = "samtools merge {}_sorted.merged.bam {}".format(name, inputs)
    """The Addgroup name function"""
    add_cmd = "gatk --java-options \"-XX:ParallelGCThreads=3 -Xmx50g\" AddOrReplaceReadGroups "
    add_cmd += "-I= {}_sorted.merged.bam ".format(name)
    add_cmd += "-O= {}_sorted.merged.addg.bam ".format(name)
    add_cmd += "-RGID={} ".format(ID)
    add_cmd += "-RGLB=lib1 "
    add_cmd += "-RGPL=BGISEQ "
    add_cmd += "-RGSM={} ".format(ID)
    add_cmd += "-RGPU=unit1 "
    ##print(add_cmd)
    """Create a .sh files with the samtools merge and the add group name functions."""
    file = open('{}_merging_addname.sh'.format(output),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem-per-cpu 50G \n')
    file.write('#SBATCH -c 3 \n')
    file.write('#SBATCH --time=24:00:00 \n')
    file.write(mg_cmd)
    file.write('\n')
    file.write('echo \" THE MERGING IS DONE ##########################################################\" \n')
    file.write(add_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}_merging_addname.out {}_merging_addname.sh".format(output, output)
    subprocess.call(sub_cmd, shell=True)

def add_only(name, ID, output):
    """No need to merged so just copy."""
    mg_cmd = "cp {}_0.sorted.bam {}_sorted.merged.bam".format(name, name)
    ##print(mg_cmd)
    """The Addgroup name function"""
    add_cmd = "gatk --java-options \"-XX:ParallelGCThreads=3 -Xmx50g\" AddOrReplaceReadGroups "
    add_cmd += "-I= {}_sorted.merged.bam ".format(name)
    add_cmd += "-O= {}_sorted.merged.addg.bam ".format(name)
    add_cmd += "-RGID={} ".format(ID)
    add_cmd += "-RGLB=lib1 "
    add_cmd += "-RGPL=BGISEQ "
    add_cmd += "-RGSM={} ".format(ID)
    add_cmd += "-RGPU=unit1 "
    ##print(add_cmd)
    """Create a .sh files with the samtools merge and the add group name functions."""
    file = open('{}_merging_addname.sh'.format(output),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem-per-cpu 50G \n')
    file.write('#SBATCH -c 3 \n')
    file.write('#SBATCH --time=12:00:00 \n')
    file.write(mg_cmd)
    file.write('\n')
    file.write('echo \" THE MERGING IS DONE ##########################################################\" \n')
    file.write(add_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}_merging_addname.out {}_merging_addname.sh".format(output, output)
    subprocess.call(sub_cmd, shell=True)


def map_sum(name, direct, output):
    """The samtools summary function."""
    sum_cmd = "samtools flagstat {} >> {}summary/{}_summary.txt".format(name, direct, output)
    ##print(sum_cmd)
    """Create a .sh files with the samtools summary function."""
    file = open('{}summary/{}_summary.sh'.format(direct, output),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 50G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=10:00:0 \n')
    file.write('echo This is information for {} > {}summary/{}_summary.txt'.format(name, direct, output))
    file.write('\n')
    file.write(sum_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}summary/{}_summary.out {}summary/{}_summary.sh".format(direct, output, direct, output)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# Merged all the lane for one sample

for name in bamfile_dir:
    print(name)
    if os.path.exists("{}/{}/bam_files/{}_sorted.merged.bam".format(path, sp, name)):
        print("The merged file for {} exists".format(name))
        if os.path.exists("{}/{}/bam_files/{}_sorted.merged.addg.bam".format(path, sp, name)):
            print("The group name has also been added --> {} IS DONE".format(name))
        else:
            print("The group name hasn't been added")
    else:
        print("The merged file for {} doesn't exist".format(name))
        if len(bamfile_dir[name][0])<2:
            print("There is only one lane for the samples --> nothing to merged, copy")
            add_only(name="{}/{}/bam_files/{}".format(path, sp, name), ID=name, output="{}/{}/bam_files/{}".format(path, sp, name))
            print("Also getting the summary of each lanes")
            for l in bamfile_dir[name][0]:
                last_path = os.path.basename(l)
                output = last_path.split('.sort')[0]
                map_sum(name=l, direct=directory, output=output)
        else:
            inputbams = ' '.join(bamfile_dir[name][0])
            if all([os.path.exists(p) for p in bamfile_dir[name][0]]):
                print("All the bam files for {} exist and will be merged".format(name))
                merge_add(name="{}/{}/bam_files/{}".format(path, sp, name), inputs=inputbams, ID=name, output="{}/{}/bam_files/{}".format(path, sp, name))
                print("Also getting the summary of each lanes")
                for l in bamfile_dir[name][0]:
                    last_path = os.path.basename(l)
                    output = last_path.split('.sort')[0]
                    map_sum(name=l, direct=directory, output=output)
            else:
                print("Oooh NoOoOo at least some directories for {} doesn't exist so {} can't be merged".format(name,name))

