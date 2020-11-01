set -e

which samtools &>/dev/null || { echo "ERROR: samtools not found!"; exit 1; }
samtools --version

which bamCoverage &>/dev/null || { echo "ERROR: deeptools not found!"; exit 1; }
bamCoverage --version

which bowtie &>/dev/null || { echo "ERROR: bowtie not found!"; exit 1; }
bowtie --version

which fastqc &>/dev/null || { echo "ERROR: fastqc not found!"; exit 1; }
fastqc --version

which macs2 &>/dev/null || { echo "ERROR: MACS2 not found!"; exit 1; }
macs2 --version

#which Rscript &>/dev/null || { echo "ERROR: R not found!"; exit 1; }
#Rscript --verion

which java &>/dev/null || { echo "ERROR: java not found!"; exit 1; }
java -version

which bedtools &>/dev/null || { echo "ERROR: bedtools not found!"; exit 1; }
bedtools --version

which SICER.sh &>/dev/null || { echo "ERROR: SICER not found!"; exit 1; }
SICER.sh --version


