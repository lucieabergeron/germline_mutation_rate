#!/bin/bash
##
PATH=$path
mkdir 1.map/submit_split
COUNT=0
echo "#!/bin/bash" > $PATH/1.map/mapping_splitted.sh
while read p; do
        COUNT=$(( $COUNT + 1 ))
        echo "#!/bin/bash" > $PATH/1.map/submit_split/submit_map_$COUNT
        echo "#SBATCH --partition normal" >> $PATH/1.map/submit_split/submit_map_$COUNT
        echo "#SBATCH --mem-per-cpu 5G" >> $PATH/1.map/submit_split/submit_map_$COUNT
        echo "#SBATCH -c 4" >> $PATH/1.map/submit_split/submit_map_$COUNT
        echo $p >> $PATH/1.map/submit_split/submit_map_$COUNT;
        ##
        echo "sbatch -o $PATH/1.map/submit_split/submit_map_$COUNT.out -e $PATH/1.map/submit_split/submit_map_$COUNT.out $PATH/1.map/submit_split/submit_map_$COUNT" >> $PATH/1.map/mapping_splitted.sh
done < $PATH/1.map/mapping.sh
echo "Split is done on $COUNT sequences"
##
