#!/bin/bash
docker run --rm --name sb -d -p 8711:8787 -v /mnt/vol1/MED_GEN_DATA:/mnt -v /mnt/vol1/student00:/home/student -t sbw-medgenetics
