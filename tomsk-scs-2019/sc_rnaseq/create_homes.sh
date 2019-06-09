#!/bin/bash
STUDENTS=26
mnt="/mnt/t/LXC_CT_data/Ubuntu1804IFMOx64/sc"
for t in `seq -w 0 $STUDENTS`; do 
    vol="/mnt/t/LXC_CT_data/Ubuntu1804IFMOx64/homes/student$t"
	echo $t $vol
    mkdir $vol
    chown 1000:1000 $vol
done

