#!/bin/bash
##
mkdir 1.map/submit_split
COUNT=0
echo "#!/bin/bash" > 1.map/mapping_splitted.sh
while read p; do
        COUNT=$(( $COUNT + 1 ))
        echo "#!/bin/bash" > 1.map/submit_split/submit_map_$COUNT
        echo "#SBATCH --partition normal" >> 1.map/submit_split/submit_map_$COUNT
        echo "#SBATCH --mem-per-cpu 5G" >> 1.map/submit_split/submit_map_$COUNT
        echo "#SBATCH -c 4" >> 1.map/submit_split/submit_map_$COUNT
        echo $p >> 1.map/submit_split/submit_map_$COUNT;
        ##
        echo "sbatch -o 1.map/submit_split/submit_map_$COUNT.out -e 1.map/submit_split/submit_map_$COUNT.out 1.map/submit_split/submit_map_$COUNT" >> 1.map/mapping_splitted.sh
done < 1.map/mapping.sh
echo "Split is done on $COUNT sequences"
##
