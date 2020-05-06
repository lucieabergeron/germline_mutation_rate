#!/bin/bash
##
## This script evaluate the insert size of the PE sequences by:
## 3. evaluating the distance between the mapping.
## 2.eva: evaluate size ##############################################################################################################################################
## Construct eva.sh with this command:
##
mkdir 2.eva
echo '#!/bin/bash' > 2.eva/eva.sh
echo '#SBATCH --partition normal' >> 2.eva/eva.sh
echo '#SBATCH --mem 16G' >> 2.eva/eva.sh
ls 1.map/map_seq | sed 's/_read.*//' | sort | uniq | while read p; do echo "evalInsertSize.pl 1.map/map_seq/"$p"_read_1.smp.bam 1.map/map_seq/"$p"_read_2.smp.bam > 2.eva/"$p".eis"; done >> 2.eva/eva.sh
##
## Then run the eva.sh:
##
sbatch 2.eva/eva.sh
##
## When this is done you can run for results
