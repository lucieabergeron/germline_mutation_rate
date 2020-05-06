y_mapping.sh
SUM=$1
OUT=$2
SEQ=$(basename $SUM | sed -e 's/\(_summary.txt\)*$//g')

PASSED=$(grep "paired in sequencing" $SUM | grep -Eo '^[0-9]+')
MAPPED=$(grep "properly paired" $SUM | grep -Eo '^[0-9]+')
echo $SEQ $PASSED . $MAPPED . >> $OUT/mapped_summary.txt
