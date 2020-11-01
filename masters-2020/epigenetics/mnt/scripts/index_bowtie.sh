#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author os@jetbrains.com.com

# Stop exec on error.
set -e

which bowtie &>/dev/null || { echo "ERROR: bowtie not found!"; exit 1; }

>&2 echo "index-bowtie $@"
if [[ $# -lt 2 ]]; then
    echo "Need 2 parameters! <GENOME> <FOLDER>"
    exit 1
fi
GENOME=$1
FOLDER=$2

cd ${FOLDER}
# Check both 32 and 64 large indexes
if ([[ ! -f "$GENOME.1.ebwt" ]] && [[ ! -f "$GENOME.1.ebwtl" ]]); then
    bowtie-build $(find . -type f -name "*.fa" | sed 's#\./##g' | paste -sd "," -) ${GENOME} 2>&1 |\
        tee ${GENOME}_bowtie_indexes.log
fi

>&2 echo "Done. index-bowtie $@"