#!/bin/bash
STUDENTS=26

for t in `seq -w 1 $STUDENTS`; do 
    docker rm sb$t.old
done
