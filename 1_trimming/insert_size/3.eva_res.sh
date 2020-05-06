#!/bin/bash
##
## This script evaluate the insert size of the PE sequences by:
## 4. formatting the results
## Got the results:
##
echo '#LANE   NUM_of_PAIRS    MEAN    MEDIAN  SD' > 2.eva/res_ins.txt
cat $(ls 2.eva/*.eis) | grep -v '^#' >> 2.eva/res_ins.txt
##
## Got ALL the results:
##
## Remove the '_read' from the file:
sed -ri 's:_read::g' 2.eva/res_ins.txt
##
## Add the column to each line and save in another file
echo $(grep '#' res_ins.txt) MAX MIN > 2.eva/res_ins_fin.txt
grep -v '#' 2.eva/res_ins.txt | while read a b c d e; do echo "$a $b $c $d $e $(echo $c + 3*$e | bc) $(echo $c - 3*$e | bc)";done >> 2.eva/res_ins_fin.txt
##
## Now the results are in res_ins_fin.txt
