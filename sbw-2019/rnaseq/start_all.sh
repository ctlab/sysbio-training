#!/bin/bash
STUDENTS=26
for t in `seq -w 0 $STUDENTS`; do 
    t10=$(echo $t | sed 's/^0//g'); 
	vol="/mnt/vol$(((t10)%4+1))/student$t"
	echo $t $vol
    docker run --name sb$t -m 16g --cpus=4 -d -p 87$t:8787 \
        -v /mnt/vol1/sbw-data/:/mnt:ro \
        -v $vol:/home/student \
        -t sbw-rnaseq 
    break
done
