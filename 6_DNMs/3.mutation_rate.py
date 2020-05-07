# -*- coding: utf-8 -*-
"""
This script calculate the mutation rate per trio and gives information on the type of mutation
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
direct_denovo="{}/{}/de_novo_mutation/".format(path, sp)
direct_handling="{}/{}/vcf_handling/".format(path, sp)
direct="{}/{}/".format(path, sp)

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

# Import nb mutation and fnr
nb_denovo = pd.read_csv('{}denovo.txt'.format(direct_denovo),sep=' ', index_col=None, header=None)
fnr = pd.read_csv('{}fnr.txt'.format(direct_denovo),sep=' ', index_col=None, header=None)

##################################################
# What you run  ##################################
##################################################
## Move the file you don't want:
list_exist=[]
for name in trio_dir:
    list_exist.append(os.path.exists("{}callability_{}.out".format(direct_denovo, name)))
list_exist.append(os.path.exists("{}fnr.txt".format(direct_denovo)))
list_exist.append(os.path.exists("{}denovo.txt".format(direct_denovo)))
if all(list_exist):
    print("Number of denovo, callability and false negative has been calculated")
    mv_log= "mv -t {}log_file/ {}callability_*.sh {}false_neg_* {}output* {}*_fnr* {}*_nb_denovo*".format(direct_denovo, direct_denovo, direct_denovo, direct_denovo, direct_denovo, direct_denovo)
    subprocess.call(mv_log, shell=True)
    print("Move log files")
    mv_vcf= "mv -t {}inter_vcf/ {}callability_filt* {}neg_rate*".format(direct_denovo, direct_denovo, direct_denovo)
    subprocess.call(mv_vcf, shell=True)
    print("Move inter vcf files")
    # Find the mutation rate per trio: ###############
    # Create a file with the callability and import:
    for name in trio_dir:
        call_cmd= "echo \"{} $(awk 'NR==4' {}callability_{}.out)\" >> {}callability.txt".format(name, direct_denovo, name, direct_denovo)
        subprocess.call(call_cmd, shell=True)
##
    call = pd.read_csv('{}callability.txt'.format(direct_denovo),sep=' ', index_col=None, header=None)
    #
    # Find the overall fnr:
    # Alpha before only AB
    fnr_all=round(sum(fnr[2])/sum(fnr[1]),5)
    # Alpha allelic balance AND site filter
    a_RP=0.002699796
    a_MQRS=0.0227818
    a_FS=0.01
    alpha_all = 1 - ((1-a_FS)*(1-a_RP)*(1-a_MQRS)*(1-fnr_all))
    print('alpha={}'.format(alpha_all))
    #
    # Find the mutation rate per trios:
    file = open('{}mutation_rate.txt'.format(direct_denovo),'w')
    for name in trio_dir:
        nb_mut=nb_denovo.loc[nb_denovo[0] ==name][1]
        C=call.loc[call[0] ==name][1]
        alpha=alpha_all
        #alpha=fnr.loc[fnr[0] ==name][3]
        mut_rate = float(nb_mut)/(2*float(C)*(1-float(alpha)))
        file.write('{} {} \n'.format(name, mut_rate))
    #
    file.close()
else:
    print("Something is wrong the denovo or callability or false negative rate hasn't been calculated")
