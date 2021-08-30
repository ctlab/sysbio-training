#!/bin/bash
STUDENTS=30
TEMPLATE="/mnt/SSDTemp/scs2021s/student-template/"
for t in `seq -w 0 $STUDENTS`; do 
    vol="/mnt/SSDTemp/scs2021s/homes/student$t"
    echo $t $vol
    mkdir $vol
    chmod 777 $vol
done
