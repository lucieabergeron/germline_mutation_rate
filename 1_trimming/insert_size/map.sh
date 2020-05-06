#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 1G
if [ $1 ]
then
    ls 0.smp/smp_seq/*.gz | while read p; do echo "bwa aln -t 4 $1 $p > 1.map/map_seq/$(basename $p '.fq.gz').sai && bwa samse $1 1.map/map_seq/$(basename $p '.fq.gz').sai $p | samtools view -Sb - > 1.map/map_seq/$(basename $p '.fq.gz').bam && rm -f 1.map/map_seq/$(basename $p '.fq.gz').sai"; done > 1.map/mapping.sh
fi
