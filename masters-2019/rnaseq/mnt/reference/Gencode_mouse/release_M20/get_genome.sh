#!/bin/bash
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz

wget  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz
gunzip gencode.vM20.annotation.gtf.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.transcripts.fa.gz

wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_Gencode_VM18.bed.gz
gunzip mm10_Gencode_VM18.bed.gz
wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_rRNA.bed.gz
gunzip mm10_rRNA.bed.gz
wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10.HouseKeepingGenes.bed.gz
gunzip mm10.HouseKeepingGenes.bed.gz

wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
gunzip refFlat.txt.gz
