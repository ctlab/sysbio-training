#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

which bamCoverage &>/dev/null || { echo "ERROR: deeptools not found!"; exit 1; }

# Load utils
source ./util.sh

>&2 echo "Batch bigwig $@"
if [[ $# -lt 1 ]]; then
    echo "Need 1 parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -vE ".tr")
do :
    NAME=${FILE%%.bam} # file name without extension
    BW=${NAME}.bw
    if [[ ! -f ${BW} ]]; then
        bamCoverage --bam ${FILE} -o ${BW} --normalizeUsing RPKM 2>&1 | tee ${NAME}_bw.log
    fi
done

check_logs

>&2 echo "Done. Batch bigwig $@"
