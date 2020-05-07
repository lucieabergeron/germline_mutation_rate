
# -*- coding: utf-8 -*-
"""
This script can be run when back has been splitted per trios and:
       - move the previous log files
       - find the mean depth save it in a file
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
direct = "{}/{}/vcf_handling/".format(path, sp)
ref="{}/{}/ref_fasta/{}.fa".format(path, sp,refGenome)
dir_tab="{}/{}/".format(path, sp)

# Dictionary:
f = open('{}/{}/pedigree.ped'.format(path, sp))
trio_dir = {}
for line in f:
    off = line.split()[1]
    fa = line.split()[2]
    mo = line.split()[3]
    name = off
    if name not in trio_dir:
        trio_dir[name] = []
    trio_dir[name].append((off, fa, mo))

# The functions:
def back_handling(direct, name):
    """Extract the depth for each samples vcf"""
    dp_1_cmd="egrep -v \"^#\" {}back_combine_genomicDBI_{}.g.vcf | cut -f 10 | sed 's/[^:]*://' | sed 's/[^:]*://' | sed 's/:.*$//' >>{}DP_{}_1.txt".format(direct, name, direct, name)
    dp_2_cmd="egrep -v \"^#\" {}back_combine_genomicDBI_{}.g.vcf | cut -f 11 | sed 's/[^:]*://' | sed 's/[^:]*://' | sed 's/:.*$//' >>{}DP_{}_2.txt".format(direct, name, direct, name)
    dp_3_cmd="egrep -v \"^#\" {}back_combine_genomicDBI_{}.g.vcf | cut -f 12 | sed 's/[^:]*://' | sed 's/[^:]*://' | sed 's/:.*$//' >>{}DP_{}_3.txt".format(direct, name, direct, name)
    """Awk to have only numbers"""
    awk1= "awk '$1 !~ /[^0-9]+/' {}DP_{}_1.txt>{}DP_{}_1_rm.txt".format(direct, name, direct, name)
    awk2= "awk '$1 !~ /[^0-9]+/' {}DP_{}_2.txt>{}DP_{}_2_rm.txt".format(direct, name, direct, name)
    awk3= "awk '$1 !~ /[^0-9]+/' {}DP_{}_3.txt>{}DP_{}_3_rm.txt".format(direct, name, direct, name)
    """Catenate the three individuals"""
    cat="cat {}DP_{}_1_rm.txt {}DP_{}_2_rm.txt {}DP_{}_3_rm.txt > {}DP_{}_tot_rm.txt".format(direct, name, direct, name, direct, name, direct, name)
    """Python script to find the mean depth per individuals and for all of them"""
    subprocess.call("cp coverage_python.py {}{}_coverage_python.py".format(direct, name), shell=True)
    subprocess.call("sed -i '1s/^/sp=\"{}\"\\n/' {}{}_coverage_python.py \n".format(sp, direct, name), shell=True)
    subprocess.call("sed -i '1s/^/path=\"{}\"\\n/' {}{}_coverage_python.py \n".format(path, direct, name), shell=True)
    subprocess.call("echo \"average(direct=direct, name='DP_{}_1_rm.txt')\" >> {}{}_coverage_python.py \n".format(name, direct, name), shell=True)
    subprocess.call("echo \"average(direct=direct, name='DP_{}_2_rm.txt')\" >> {}{}_coverage_python.py \n".format(name, direct, name), shell=True)
    subprocess.call("echo \"average(direct=direct, name='DP_{}_3_rm.txt')\" >> {}{}_coverage_python.py \n".format(name, direct, name), shell=True)
    subprocess.call("echo \"average(direct=direct, name='DP_{}_tot_rm.txt')\" >> {}{}_coverage_python.py \n".format(name, direct, name), shell=True)
    """Create a .sh files with all the functions."""
    file = open('{}back_handling_{}.sh'.format(direct, name),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --mem 150G \n')
    file.write('#SBATCH -c 20 \n')
    file.write('#SBATCH --time=10:00:00 \n')
    file.write('## Extract depth per samples \n')
    file.write(dp_1_cmd)
    file.write('\n')
    file.write(dp_2_cmd)
    file.write('\n')
    file.write(dp_3_cmd)
    file.write('\n')
    file.write('## Remove non number with awk \n')
    file.write(awk1)
    file.write('\n')
    file.write(awk2)
    file.write('\n')
    file.write(awk3)
    file.write('\n')
    file.write('## Catenate the three individuals \n')
    file.write(cat)
    file.write('\n')
    file.write('## Calculate the median \n')
    file.write('python {}{}_coverage_python.py'.format(direct, name))
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}depth_{}.out {}back_handling_{}.sh".format(direct, name, direct, name)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################
# Filter per trio files:
list_exist=[]
for name in trio_dir:
    list_exist.append(os.path.exists("{}back_combine_genomicDBI_{}.g.vcf.idx".format(direct, name)))
if all(list_exist):
    print("Back combine VCF files for each trio exists --> find the depth")
    mv_log= "mv -t {}log_file/ {}back*.out {}back*.sh".format(direct, direct, direct)
    subprocess.call(mv_log, shell=True)
    print("Move the log files")
    for name in trio_dir:
        back_handling(direct=direct, name=name)
        print("Depth of back combine VCF file for {} --> depth_{}.out".format(name, name))
else:
    print("Some back combine VCF files per trios doesn't exists")
