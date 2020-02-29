#!/bin/bash
t=$1
tag=sbw-medgenetics
if [ -z "$t" ]; then
    echo give a number
fi

t10=$(echo $t | sed 's/^0//g'); 
vol="/mnt/vol$(((t10)%4+1))/student$t"
echo $t $vol
docker run --name sb$t -m 16g --cpus=4 -d -p 87$t:8787 \
    -v /mnt/vol1/sbw-data/:/mnt:ro \
    -v $vol:/home/student \
    -t "$tag"
