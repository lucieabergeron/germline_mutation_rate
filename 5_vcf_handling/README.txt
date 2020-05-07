Set up the right variable in variable.py:

Test relatedness, create a vcf per trio, filter the snps (site filter), select mendelian violations and extract the average depth (from vcf):
    python 0.relate_split_trios.py --> check the relatedness between individuals and split genotypes/back_combine per trios
    python 1.geno_splitted_handling.py --> filter the sites, find the mendelian violation and convert to table
    python 2.back_splitted_handling.py --> depth
    python 3.final_check.py --> check if everything is here and move the files
