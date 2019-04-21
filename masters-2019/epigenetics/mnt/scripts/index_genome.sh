#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

>&2 echo "index-genome $@"
if [[ $# -lt 2 ]]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2


if [[ ! -d "$FOLDER" ]]; then
    echo "Created indexes folder: ${FOLDER}"
    mkdir -p ${FOLDER}
fi
cd ${FOLDER}

echo "Check chrom.sizes file"
if [[ ! -f "${GENOME}.chrom.sizes" ]]; then

    echo "Download only chromosomes sequences"
    rsync -avzP --exclude="*.txt" rsync://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/chromosomes/ .
    gunzip -f *.fa.gz

    echo "Downloading chrom.sizes file"
    wget -nc http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.chrom.sizes

    chmod a+r *
fi

>&2 echo "Done. index-genome $@"
