#!/bin/bash
t=$1
tag=sb-scrnaseq
if [ -z "$t" ]; then
    echo give a number
    exit 0
fi

vol="/mnt/t/LXC_CT_data/Ubuntu1804IFMOx64/homes/student$t"
mnt="/mnt/t/LXC_CT_data/Ubuntu1804IFMOx64/sc"
echo $t $vol
docker run --name sb$t -m 16g -d -p 87$t:8787 \
    -v $mnt:/mnt \
    -v $vol:/home/student \
    -t "$tag"

