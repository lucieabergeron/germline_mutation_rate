Set up the right variable in variable.py

Post mapping processing and quality checks:
    python 0.uniq_rmdup.py --> keep only reads mapping a unique region of the genome, remove duplicates and indexes the files (require samtools version 1.4.1)
    python 1.coverage.py --> give the coverage and number of reads (require samtools version 1.2)
    python 2.final_check.py --> see if all files are here and remove intermediate files
