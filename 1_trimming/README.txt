Set up the right variable in variable.py:

evaluate insert size in insert_size/
    0.smp.sh --> sub samples the sequences
    1.map.sh --> SE mapping (require bwa version 0.7.15 and samtools version 1.2)
    2.eva.sh --> evaluate insert size (require samtools version 0.1.19 and the directory to it in evalInsertSize.pl)
    3.eva_res.sh --> format the results

Trimmed the sequences in filtering/
    0.filtering.sh --> filters the sequences (require soapnuke version: 1.5.6)
    1.summary_filtering.sh --> to get the summary
