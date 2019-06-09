#!/bin/bash
STUDENTS=26
mnt="/mnt/t/LXC_CT_data/Ubuntu1804IFMOx64/sc"
tag=sb-scrnaseq

for t in `seq -w 1 $STUDENTS`; do 
    vol="/mnt/t/LXC_CT_data/Ubuntu1804IFMOx64/homes/student$t"
	echo $t $vol
    docker run --name sb$t -m 8g --cpus=4 -d -p 87$t:8787 \
        -v $mnt/:/mnt:ro \
        -v $vol:/home/student \
        -t $tag
done
