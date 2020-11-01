Epigenetics practice
====================

The materials from the Epigenetics module on November 2-7, 2020. 

By Oleg Shpynov(os@jetbrains.com) and Roman Chernyatchik(roman.chernyatchik@jetbrains.com).



Build Docker Image
------------------

```bash
docker build -t jbrepi .
```



Launch instance
----------------------

Append `/bin/bash` and check that all the tools are installed correctly with `check.sh` script. 
```
docker run --rm --name student1 -m 16g --cpus=2 -p 8787:8787 -t jbrepi
```