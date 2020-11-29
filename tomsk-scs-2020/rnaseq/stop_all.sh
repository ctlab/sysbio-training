#!/bin/bash
STUDENTS=26
for t in `seq -w 0 $STUDENTS`; do 
    docker stop sb$t
    docker rm sb$t
done
