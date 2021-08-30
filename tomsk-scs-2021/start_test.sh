#!/bin/bash
t=00
shared="/mnt/SSDTemp/scs2021s/shared"
vol="/mnt/SSDTemp/scs2021s/homes/student$t"
docker run --name sb$t -m 4g --cpus=4 -d -p 87$t:8787 \
    -v $vol:/home/student \
    -v $shared:/home/student/shared:ro \
    -t ctlab.registry.jetbrains.space/p/bitrain/docker/sysbio-2021
