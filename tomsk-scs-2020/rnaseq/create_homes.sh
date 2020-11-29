#!/bin/bash
STUDENTS=26
for t in `seq -w 0 $STUDENTS`; do 
    t10=$(echo $t | sed 's/^0//g'); 
	vol="/mnt/vol$(((t10)%4+1))/student$t"
	echo $t $vol
    mkdir $vol
    chown 1000:1000 $vol
done

