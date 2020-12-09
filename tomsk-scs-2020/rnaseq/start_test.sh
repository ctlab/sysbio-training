#!/bin/bash
docker run --rm --name sb -d -p 8700:8787 -v /mnt/scs2020/students/student00:/home/student -v /mnt/scs2020/shared:/home/student/shared:ro -t tomsk-scs-2020-rnaseq 

