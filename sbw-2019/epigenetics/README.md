SBW-2019 epigenetics
====================

The materials from the SBW-2019 Epigenetics module on May 22, 2019. 

By Oleg Shpynov(oleg.shpynov@jetbrains.com) and Roman Chernyatchik(roman.chernyatchik@jetbrains.com).


Build Docker Image
------------------

```bash
docker build -t sbw-epi .
```

Batch launch for students at IFMO infrastucture
-----------------------------------------------
Important: student home folders are split across different HD volumnes `/mnt/vol1` - `/mnt/vol4`.

```bash
SS=1
SE=26
for t in `seq -w $SS $SE`; do 
    t10=$(echo $t | sed 's/^0//g');  
	vol=/mnt/vol$(((t10)%4+1))/student$t
	echo $t $vol
    docker run --name sb$t -m 32g --cpus=6 -d -p 87$t:8787 \
        -v /scratch/oshpynov/mnt:/mnt:ro \
        -v /mnt/vol$(((t10)%4+1))/student$t/student$t:/home/student \
        -t sbw-epi; 
done
```

Open localhost:871 and use credential student:sysbiopass to login. 

Stop containers
---------------

```bash
#!/bin/bash
STUDENTS=26
for t in `seq -w 0 $STUDENTS`; do 
    docker stop sb$t
    docker rm sb$t
done
```

Cleanup Docker
--------------

```bash
docker rmi $(docker images --filter "dangling=true" -q --no-trunc)
```