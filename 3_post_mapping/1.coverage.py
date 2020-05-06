# -*- coding: utf-8 -*-
"""
This script gather information about the coverage of the bam files.
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
in_dir = "{}/{}/bam_files/".format(path, sp)
out_dir = "{}/{}/bam_files/coverage/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}_sorted.merged.addg.bam".format(path, sp, name)


# The function:
def coverage(in_dir, seq, out_dir):
    """The samtools function"""
    cov_cmd = "samtools depth {}{}_sorted.merged.addg.uniq.rmdup.bam 1> {}{}_coverage.txt".format(in_dir, seq, out_dir, seq)
    print("\t samtools depth for coverage --> {}{}_coverage.txt".format(out_dir, seq))
    """Summarize coverage"""
    sum_cmd = "python {}{}_coverage_python.py".format(out_dir, seq)
    print("\t Summarize the depth with a python script --> {}summary_coverage.txt".format(out_dir))
    """Prepare"""
    subprocess.call("cp coverage_python.py {}{}_coverage_python.py".format(out_dir, name), shell=True)
    subprocess.call("echo \"average_column(csv='{}{}_coverage.txt', name='{}')\" >> {}{}_coverage_python.py \n".format(out_dir, name, name, out_dir, name), shell=True)
    subprocess.call("sed -i '1s/^/sp=\"{}\"\\n/' {}{}_coverage_python.py \n".format(sp, out_dir, name), shell=True)
    """Create a .sh files."""
    file = open('{}{}_coverage.sh'.format(out_dir, seq),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 30G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=06:00:0 \n')
    file.write(cov_cmd)
    file.write('\n')
    file.write(sum_cmd)
    file.write('\n')
    file.close()
    """Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}{}_coverage.out {}{}_coverage.sh".format(out_dir, seq, out_dir, seq)
    subprocess.call(sub_cmd, shell=True)

# The function for getting the number of reads:
def nbread_sum(seq, direct):
    """The samtools summary function."""
    nb_merged_cmd = "samtools view {}{}_sorted.merged.addg.bam | wc -l".format(direct, seq)
    nb_uniq_cmd = "samtools view {}{}_sorted.merged.addg.uniq.bam | wc -l".format(direct, seq)
    nb_rmdup_cmd = "samtools view {}{}_sorted.merged.addg.uniq.rmdup.bam | wc -l".format(direct, seq)
    print("\t Run number of reads in summary")
    """Create a .sh files with the samtools summary function."""
    file = open('{}summary/{}_nb_reads.sh'.format(direct, seq),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 50G \n')
    file.write('#SBATCH -c 1 \n')
    file.write('#SBATCH --time=10:00:0 \n')
    file.write('SEQ={} \n'.format(seq))
    file.write('NB_MERGED=$({}) \n'.format(nb_merged_cmd))
    file.write('NB_UNIQ=$({}) \n'.format(nb_uniq_cmd))
    file.write('NB_RMDUP=$({}) \n'.format(nb_rmdup_cmd))
    file.write('echo $SEQ $NB_MERGED . $NB_UNIQ $NB_RMDUP >> {}summary/nb_reads.txt \n'.format(direct))
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}summary/{}_nb_reads.out {}summary/{}_nb_reads.sh".format(direct, seq, direct, seq)
    subprocess.call(sub_cmd, shell=True)

# The function for getting the name of chromosomes:
def chr_name(seq, direct):
    """The samtools summary function."""
    chr_cmd = "samtools idxstats {}{}_sorted.merged.addg.uniq.rmdup.bam > {}summary/chrom_name.txt".format(direct, seq, direct)
    """Create a .sh files with the chrom name function."""
    file = open('{}summary/chrom_name.sh'.format(direct),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal \n')
    file.write('#SBATCH --mem 64G \n')
    file.write('#SBATCH -c 16 \n')
    file.write('#SBATCH --time=10:00:0 \n')
    file.write(chr_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}summary/chrom_name.out {}summary/chrom_name.sh".format(direct, direct)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################
for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}_sorted.merged.addg.uniq.rmdup.bam".format(in_dir, name)):
        print("\t The uniq_rmdup file for {} exists".format(name))
        # Count the number of reads:
        nbread_sum(seq=name, direct=in_dir)
    if os.path.exists("{}{}_coverage.txt".format(out_dir, name)):
        print("\t The coverage of {} has already been calculated --> NOT CALCULATED".format(name))
    else:
        print("\t The coverage of {} has NOT been calculated".format(name))
        coverage(in_dir=in_dir, seq=name, out_dir=out_dir)
        # Prepare the python script
##        subprocess.call("cp coverage_python.py {}{}_coverage_python.py".format(out_dir, name), shell=True)
##        subprocess.call("echo \"average_column(csv='{}{}_coverage.txt', name='{}')\" >> {}{}_coverage_python.py \n".format(out_dir, name, name, out_dir, name), shell=True)
##        subprocess.call("sed -i '1s/^/sp=\"{}\"\\n/' {}{}_coverage_python.py \n".format(sp, out_dir, name))

chr_name(seq=name, direct=in_dir)
print("\n Look at the name of chromosome only for {} sequence".format(name))
