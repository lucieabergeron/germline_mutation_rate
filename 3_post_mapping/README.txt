Set up the right variable in variable.py

Post mapping processing and quality checks:
    python 0.uniq_rmdup.sh --> keep only reads mapping a unique region of the genome, remove duplicates and indexes the files (require samtools version 1.4.1)
    python 1.coverage.sh --> give the coverage and number of reads (require samtools version 1.2)
    python 2.final_check.sh --> see if all files are here and remove intermediate files
