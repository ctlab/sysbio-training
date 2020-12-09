#!/bin/bash
STUDENTS=26
for t in `seq -w 0 $STUDENTS`; do 
    t10=$(echo $t | sed 's/^0//g'); 
	vol="/mnt/scs2020/students/student$t"
	echo $t $vol
    rsync -raP /mnt/scs2020/students/template/ $vol/
done

