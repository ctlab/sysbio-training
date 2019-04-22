#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

which java &>/dev/null || { echo "ERROR: java not found!"; exit 1; }

# Load utils
source $(dirname $0)/util.sh

>&2 echo "Batch remove_duplicates $@"
if [[ $# -lt 1 ]]; then
    echo "Need 2 parameters! <PICARD_TOOLS_JAR> <WORK_DIR>"
    exit 1
fi

# Check Picard tools
PICARD_TOOLS_JAR=$1
if [[ ! -f "${PICARD_TOOLS_JAR}" ]]; then
    echo "Picard tools not found! Download Picard: <http://broadinstitute.github.io/picard/>"; exit 1;
fi

WORK_DIR=$2
cd ${WORK_DIR}

for FILE in $(find . -name '*.bam' | grep -v _unique | sed 's#\./##g')
do :
    NAME=${FILE%%.bam}
    UNIQUE_BAM=${NAME}_unique.bam
    METRICS=${NAME}_metrics.txt

    if [[ ! -f ${UNIQUE_BAM} ]]; then
        java -Xmx8G -jar ${PICARD_TOOLS_JAR} \
            MarkDuplicates REMOVE_DUPLICATES=true INPUT=${FILE} OUTPUT=${UNIQUE_BAM} M=${METRICS} 2>&1 |\
            tee ${NAME}_unique.log
    fi
done

check_logs

>&2 echo "Done. Batch remove_duplicates $@"