# -*- coding: utf-8 -*-
"""
This script can be run when genotype has been splitted per trios and:
       - move the previous log files
       - filter the sites
       - detect MV
       - convert to table
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
    file = open('{}pedigree_{}.ped'.format(direct, name),'w')
    file.write(line)
    file.close()

# The functions:
def geno_handling(ref, geno_file, direct, output_snp, output_filt, pedigree, output_mv, output_table, name):
    """Select SNPs on the genotype vcf"""
    snp_cmd = "gatk --java-options \"-XX:ParallelGCThreads=16 -Xmx60g \" SelectVariants "
    snp_cmd += "-R={} ".format(ref)
    snp_cmd += "-V={} ".format(geno_file)
    snp_cmd += "--select-type-to-include SNP "
    snp_cmd += "-O={}{}".format(direct, output_snp)
    """Site filter SNPs on the genotype vcf"""
    sfilt_cmd = "gatk --java-options \"-XX:ParallelGCThreads=16 -Xmx60g \" VariantFiltration "
    sfilt_cmd += "-R={} ".format(ref)
    sfilt_cmd += "-V={}{} ".format(direct, output_snp)
    sfilt_cmd += "--filter-expression \"QD < 2.0 || FS > 20.0 || MQ < 40.0 || MQRankSum < -2.0 || MQRankSum > 4.0 || ReadPosRankSum < -3.0 || ReadPosRankSum > 3.0 || SOR > 3.0\"  "
    sfilt_cmd += "--filter-name \"my_snp_filter\"  "
    sfilt_cmd += "-O={}{}".format(direct, output_filt)
    """Select mendelian violations"""
    mv_cmd = "gatk --java-options \"-XX:ParallelGCThreads=16 -Xmx60g \" SelectVariants "
    mv_cmd += "-R={} ".format(ref)
    mv_cmd += "-V={}{} ".format(direct, output_filt)
    mv_cmd += "-ped {} ".format(pedigree)
    mv_cmd += "--mendelian-violation "
    mv_cmd += "-O={}{}".format(direct, output_mv)
    """Variants to table function"""
    vtt_cmd = "gatk --java-options \"-XX:ParallelGCThreads=16 -Xmx60g \" VariantsToTable "
    vtt_cmd += "-V={}{} ".format(direct, output_mv)
    vtt_cmd += "-F CHROM -F POS -F TYPE -F REF -F ALT -F FILTER -GF GT -GF AD -GF DP -GF GQ -GF PL -GF SAC "
    vtt_cmd += "-O={}{} ".format(direct, output_table)
    vtt_cmd += "--show-filtered"
    """Create a .sh files with the filter functions."""
    file = open('{}geno_handling_{}.sh'.format(direct, name),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --partition express,normal \n')
    file.write('#SBATCH --mem 64G \n')
    file.write('#SBATCH -c 16 \n')
    file.write('#SBATCH --time=00:50:00 \n')
    file.write('## Select the SNPs \n')
    file.write(snp_cmd)
    file.write('\n')
    file.write('## Filter the SNPs \n')
    file.write(sfilt_cmd)
    file.write('\n')
    file.write('## Detect Mendelien violations \n')
    file.write(mv_cmd)
    file.write('\n')
    file.write('## Convert to table \n')
    file.write(vtt_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}geno_handling_{}.out {}geno_handling_{}.sh".format(direct, name, direct, name)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# Filter per trio files:
list_exist=[]
for name in trio_dir:
    list_exist.append(os.path.exists("{}genotype_genomicDBI_{}.g.vcf.idx".format(direct, name)))
if all(list_exist):
    print("Genotype VCF files for each trio exists --> genotype splitted handling")
    vcf_tables_directories = open("{}vcf_tables_directories.txt".format(dir_tab), "w")
    vcf_dir=[]
    for name in trio_dir:
      ## move the log files:
        mv_log= "mv -t {}log_file/ {}genotype_genomicDBI_{}.out {}genotype_genomicDBI_{}.sh".format(direct, direct, name, direct, name)
        subprocess.call(mv_log, shell=True)
        geno_handling(ref=ref, geno_file="{}genotype_genomicDBI_{}.g.vcf".format(direct, name), direct=direct, output_snp="genotype_genomicDBI_{}_snp.g.vcf".format(name), output_filt="genotype_genomicDBI_{}_snp_filt.g.vcf".format(name), pedigree="{}pedigree_{}.ped".format(direct, name), output_mv="genotype_genomicDBI_{}_snp_filt_MV.g.vcf".format(name), output_table="trio_{}_vcf_table.out".format(name), name=name)
        print("Genotype VCF file filtered, MV and convert to table for trio {} --> trio_{}_vcf_table.out".format(name, name))
        vcf_dir = [name, trio_dir[name][0][1],trio_dir[name][0][2]]
        vcf_direct = "{}trio_{}_vcf_table.out".format(direct, name)
        vcf_dir.append(vcf_direct)
        vcf_dir = ' '.join(vcf_dir)
        vcf_tables_directories.write(vcf_dir + "\n")
    vcf_tables_directories.close()
    print("All tables directories are in {}vcf_tables_directories.txt".format(dir_tab))
    print("Move the log files")
    print("{} files filtered, MV, convert to table".format(len(trio_dir)))
else:
    print("Some VCF files per trios doesn't exists")
