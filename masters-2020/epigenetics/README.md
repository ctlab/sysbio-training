Epigenetics practice
====================

The materials from the Epigenetics module on November 2-7, 2020. 

By Oleg Shpynov(os@jetbrains.com) and Roman Chernyatchik(roman.chernyatchik@jetbrains.com).



Build Docker Image
------------------

```bash
docker build -t biolabs/itmo2020 .
```



Launch instance
----------------------

# Use --restart to ensure that docker will restart the container after machine is on.
```
docker run --name student --restart=always -m 16g --cpus=2 -v /mnt:/mnt:ro  -p 8787:8787 -d -t biolabs/itmo2020
```

Additional scripts
------------------
* `prepare.sh` - prepare ChIP-seq data
* `check.sh` - check ChIP-seq tools
* `check_R.sh` - check ChIP-seq R packages