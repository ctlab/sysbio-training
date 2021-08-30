#!/bin/bash
STUDENTS=30
TEMPLATE="/mnt/SSDTemp/scs2021s/students-template/"
for t in `seq -w 0 $STUDENTS`; do 
    vol="/mnt/SSDTemp/scs2021s/homes/student$t/"
    echo $t $vol
    rsync -r $TEMPLATE $vol
    chmod -R a+rwX -R $vol
done

