#!/bin/bash
docker run --rm --name sbrw -d -p 8699:8787 -v /mnt/scs2020/shared:/mnt:rw -v /mnt/scs2020/students/student00:/home/student -t tomsk-scs-2020-rnaseq 

