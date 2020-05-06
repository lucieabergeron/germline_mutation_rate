# -*- coding: utf-8 -*-
"""
This script check if the mapping has created the files,
merged the lane from the same library,
and add a group name to the merged mapped bamfile.
"""
##################################################
# What you need ##################################
##################################################
# Packages:
import os
import subprocess
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


##################################################
# What you run  ##################################
##################################################
for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}_sorted.merged.addg.bam".format(directory, name)):
        print("The final file for {} exists ---> DONE".format(name))
        # Move the merged file
        mv_mer_cmd = "mv {}{}_sorted.merged.bam {}inter_bam/".format(directory, name, directory)
        subprocess.call(mv_mer_cmd, shell=True)
        # Move the log files
        mv_log1 = "mv -t {}map.log/ {}{}_?.sorted.out {}{}_?.sorted.sh".format(directory, directory, name, directory, name)
        subprocess.call(mv_log1, shell=True)
        mv_log2 = "mv -t {}merge.log/ {}{}_merging_addname.out {}{}_merging_addname.sh".format(directory, directory, name, directory, name)
        subprocess.call(mv_log2, shell=True)
        for i in range(len(bamfile_dir[name][0])):
            files = bamfile_dir[name][0][i]
            # Move the singular bam files
            mv_cmd = "mv {} {}inter_bam/".format(files, directory)
            subprocess.call(mv_cmd, shell=True)
            # Get the summary in good format
            last_path = os.path.basename(files)
            output = last_path.split('.sort')[0]
            sum_cmd = "bash {}/summary_mapping.sh {}summary/{}_summary.txt {}summary".format(path, directory, output, directory)
            subprocess.call(sum_cmd, shell=True)
    else:
        print("The final file for {} doesn't exist".format(name))
        print("SOMETHING IS WRONG!!!!")


