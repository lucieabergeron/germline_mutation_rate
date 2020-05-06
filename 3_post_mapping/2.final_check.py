# -*- coding: utf-8 -*-
"""
This script checked if:
 the uniq_rmdup exists and is indexed
 the coverage on this files has been done
 the recalibration files exists and is indexed
 move all non utils files (especially from the recal and index)
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
bam_dir = "{}/{}/bam_files/".format(path, sp)
cov_dir = "{}/{}/bam_files/coverage/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}".format(path, sp, name)



##################################################
# What you run  ##################################
##################################################

for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}_sorted.merged.addg.uniq.rmdup.bam".format(bam_dir, name)):
        print("\t Rmdup and uniq exists for {} ".format(name))
        if os.path.exists("{}{}_sorted.merged.addg.uniq.rmdup.bam.bai".format(bam_dir, name)):
            print("\t and has been indexed ".format(name))
            mv_idx = "mv {}{}_uniq_rmdup_idx* {}uniq_rmdup.log/".format(bam_dir, name, bam_dir)
            subprocess.call(mv_idx, shell=True)
            mv_bam1 = "mv {}{}_sorted.merged.addg.bam {}inter_bam/".format(bam_dir, name, bam_dir)
            subprocess.call(mv_bam1, shell=True)
            mv_bam2 = "mv {}{}_sorted.merged.addg.uniq.bam {}inter_bam/".format(bam_dir, name, bam_dir)
            subprocess.call(mv_bam2, shell=True)
            print("\t Move the log files and the intermediate files")
        else:
            print("\t the rmdup uniq has NOT BEEN INDEXED ".format(name))
    else:
        print("\t NO RMDUP UNIQ FOR {}".format(name))
    if os.path.exists("{}{}_sorted.merged.addg.uniq.rmdup.recal.bam".format(bam_dir, name)):
        print("\t Recalibrated bam exists for {} ".format(name))
        if os.path.exists("{}{}_sorted.merged.addg.uniq.rmdup.recal.bai".format(bam_dir, name)):
            print("\t and has been indexed ".format(name))
            mv_recal = "mv {}{}_recal* {}recal.log/".format(bam_dir, name, bam_dir)
            subprocess.call(mv_recal, shell=True)
            print("\t Move the recalibration log files")
        else:
            print("\t the recalibrated has NOT BEEN INDEXED ".format(name))
    else:
        print("\t NO RECALIBRATED FOR {}".format(name))
    if os.path.exists("{}{}_coverage.txt".format(cov_dir, name)) and os.path.exists("{}{}_coverage.sh".format(cov_dir, name)):
        print("\t Coverage has probably been calculated")
        mv_cov = "mv {}{}_coverage* {}cov.log".format(cov_dir, name, cov_dir)
        subprocess.call(mv_cov, shell=True)
        print("\t Move the coverage log files")
    else:
        print("\t Coverage has probably NOT BEEN calculated")
