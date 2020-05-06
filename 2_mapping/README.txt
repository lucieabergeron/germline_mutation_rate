Set up the right variable in variable.py:

Mapping, mergind and summary of the sequences.
   python 0.mapping.py --> map the trimmed sequences to the reference genome (require bwa 0.7.15 and samtools 1.2)
   python 1.merging.py --> merged all the lane from an individual together and add group names (require samtools 1.2 and gatk 4.0.7.0)
   python 2.final_check.py --> check all the files are there and summarized the mapping in $path/bam_files/summary/mapped_summary.txt
