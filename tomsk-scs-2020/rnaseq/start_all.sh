#!/bin/bash
STUDENTS=26
for t in `seq -w 0 $STUDENTS`; do 
    t10=$(echo $t | sed 's/^0//g'); 
	vol="/mnt/scs2020/students/student$t"
	echo $t $vol
    docker run --name sb$t -m 6g --cpus=4 -d -p 87$t:8787 \
        -v $vol:/home/student \
        -v /mnt/scs2020/shared/:/home/student/shared:ro \
        -t tomsk-scs-2020-rnaseq
done
