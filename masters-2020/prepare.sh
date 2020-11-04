# Download workshop data & indexes
sudo wget https://artyomovlab.wustl.edu/publications/supp_materials/4Oleg/2020_ITMO_epigenetics_practice/mnt.tar.gz
  tar xvf /mnt.tar.gz && mv mnt/* /mnt/ && rm /mnt.tar.gz

# Download tools
sudo wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar -O /mnt/picard.jar
sudo wget https://download.jetbrains.com/biolabs/span/span-0.13.5244.jar -O /mnt/span-0.13.5244.jar
sudo wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes -O /mnt/hg19.chrom.sizes
sudo wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gtf.gz -O /mnt/gencode.v30lift37.annotation.gtf.gz
sudo wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz -O /mnt/hg19.fa.gz

# Create references to data
mkdir ~/data && for F in $(find /mnt/data/ -type f); do ln -s $F ~/data; done

