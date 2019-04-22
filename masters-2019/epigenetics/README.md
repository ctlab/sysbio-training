Docker Image for Epigenetics part
=================================

This is just an image `ubuntu:latest` - Ubuntu LTS with all the environment setup.

Build
-----
```bash
docker build -t sbepi .
```


Launch single instance
-----------------------
```bash
DATA_FOLDER=/mnt
STUDENT_FOLDER=/tmp/student1
docker run -t sbepi --name sbepi1 -m 8g --cpus=2 -d -p 871:8787 \
    -v ${DATA_FOLDER}:/mnt:ro \
    -v ${STUDENT_FOLDER}:/home/student     
```


Batch launch for students
-------------------

```bash
STUDENTS=16
DATA_FOLDER=/mnt
STUDENT_FOLDER_PREFIX=/mnt/vol
for t in $(seq -w 1 ${STUDENTS}); do  
    docker run --name sbepi$t -m 32g --cpus=6 -d -p 87$t:8787 \
    -v ${DATA_FOLDER}:/mnt:ro -v ${STUDENT_FOLDER_PREFIX}"$(((t-1)/4+1))"/student$t:/home/student \
    -t sbepi
done
```