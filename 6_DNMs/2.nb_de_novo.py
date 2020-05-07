# -*- coding: utf-8 -*-
"""
This script write the R script to find the number of de novo mutation:
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
direct_denovo="{}/{}/de_novo_mutation/".format(path, sp)
direct_handling="{}/{}/vcf_handling/".format(path, sp)
direct = "{}/{}/".format(path, sp)


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
def nb_dn(direct, direct_denovo, direct_handling, name, GQ_lim, DP_min, DP_max, AB_max, AB_min):
    """Write the start of the Rscript"""
    file = open('{}{}_nb_denovo.r'.format(direct_denovo, name),'w')
    file.write('# Sample: \n')
    file.write('name=\"{}\" \n'.format(name))
    file.write('\n')
    file.write('# Directories: \n')
    file.write('direct_denovo=\"{}\" \n'.format(direct_denovo))
    file.write('direct_handling=\"{}\" \n'.format(direct_handling))
    file.write('direct=\"{}\" \n'.format(direct))
    file.write('\n')
    file.write('# Import mean depth: \n')
    file.write('dp_table <- read.csv(paste0(direct, \"depth.txt\"), sep =\" \", header=FALSE) \n')
    file.write('mean_dp=dp_table[which(dp_table[,1]==name),2] \n')
    file.write('\n')
    file.write('# Filters: \n')
    file.write('GQ_lim={} \n'.format(GQ_lim))
    file.write('DP_min={} \n'.format(DP_min))
    file.write('DP_max={} \n'.format(DP_max))
    file.write('AB_max={} \n'.format(AB_max))
    file.write('AB_min={} \n'.format(AB_min))
    file.write('\n')
    file.close()
    """Catenate with the rest"""
    cat_cmd = "cat rscript_nb_denovo.R >>{}{}_nb_denovo.r".format(direct_denovo, name)
    """Run r script"""
    r_cmd = "Rscript {}{}_nb_denovo.r \n".format(direct_denovo, name)
    """Create a .sh files with the filter functions."""
    file = open('{}{}_nb_denovo.sh'.format(direct_denovo, name),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition normal,short \n')
    file.write('#SBATCH --mem 10G \n')
    file.write('#SBATCH -c 6 \n')
    file.write('#SBATCH --time=10:00:00 \n')
    file.write('## Catenate the Rscript for samples and remaining part\n')
    file.write(cat_cmd)
    file.write('\n')
    file.write('## Run the Rscript: \n')
    file.write(r_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}{}_nb_denovo.out {}{}_nb_denovo.sh".format(direct_denovo, name, direct_denovo, name)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# Find the number of mutation per trios:
for name in trio_dir:
    nb_dn(direct=direct, direct_denovo=direct_denovo, direct_handling=direct_handling, name=name, GQ_lim=GQ_lim, DP_min=DP_min, DP_max=DP_max, AB_max=AB_max, AB_min=AB_min)
    print("Finding the number of mutation for trio {}".format(name))

