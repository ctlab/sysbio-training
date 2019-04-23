Epigenetics practice
====================

The materials from the Epigenetics module on April 23--25, 2019. 

By Oleg Shpynov(oleg.shpynov@jetbrains.com) and Roman Chernyatchik(roman.chernyatchik@jetbrains.com).



Build Docker Image
------------------

```bash
docker build -t sbepi .
```


Launch single instance
----------------------
```bash
DATA_FOLDER=/mnt
STUDENT_FOLDER=/tmp/student1
docker run --name sbepi1 -m 32g --cpus=6 -d -p 871:8787 \
    -v ${DATA_FOLDER}:/mnt:ro \
    -v ${STUDENT_FOLDER}:/home/student \
     -t sbepi    
```

Open localhost:871 and use credential student:sysbiopass to login. 

Batch launch for students at IFMO infrastucture
-----------------------------------------------
Important: student home folders are split across different HD volumnes `/mnt/vol1` - `/mnt/vol4`.

```bash
STUDENTS=16
for t in `seq -w 1 $STUDENTS`; do 
    t10=$(echo $t | sed 's/^0//g'); 
    docker run --name sb$t -m 32g --cpus=6 -d -p 87$t:8787 \
        -v /scratch/oshpynov/mnt:/mnt:ro \ 
        -v /mnt/vol"$(((t10-1)/4+1))"/student$t:/home/student \
        -t sbepi; 
done
```