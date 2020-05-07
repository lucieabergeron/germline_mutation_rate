Set up the right variable in variable.py

Calling variant for each individual per chromosomes, combine all individuals, joint genotype and merge the chromosomes files:
    python 0.call_res_div.sh --> call variable in bp resolution mode for each samples per chromosomes
    python 1.combine_gvcf_div.sh --> combine all the samples together per chromosomes (the output GenomicDBI database is wrote on the $SLURM_JOBID node and then copy into the current directory, depending on the machine use this might be changed to write it directely on the current directory)
    python 2.genotype_gvcf.sh --> genotype per chromosomes and back combine per chromosome (GenomicDBI database is first copy to the $SLURM_JOBID node, this might be changed)
    python 3.gather_gvcf.sh --> to have a unique genotype and back combine for all chromosomes
    python 4.final_check.py --> check all the files are present and move intermediate files
