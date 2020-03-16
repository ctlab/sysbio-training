#!/bin/bash

oldId=`docker ps | grep gene-expression | awk '{print $1;}'`
docker stop $oldId
docker pull ctlab.registry.jetbrains.space/bioinfdocker/masters-gene-expression
docker run -d --restart unless-stopped -p 8787:8787 -v /data/expression:/mnt/data:ro ctlab.registry.jetbrains.space/bioinfdocker/masters-gene-expression
