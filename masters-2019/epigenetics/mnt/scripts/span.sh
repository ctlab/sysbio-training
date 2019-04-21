#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

which java &>/dev/null || { echo "ERROR: java not found!"; exit 1; }

# Load utils
source ./util.sh

>&2 echo "Batch SPAN $@"
if [[ $# -lt 5 ]]; then
    echo "Need >= 5 parameters! <SPAN_JAR_PATH> <WORK_DIR> <GENOME> <CHROM_SIZES> <BIN> [<FDR> <GAP>]"
    exit 1
fi

SPAN_JAR_PATH=$1
if [[ ! -f "${SPAN_JAR_PATH}" ]]; then
    >&2 echo "SPAN not found! Download SPAN: <https://research.jetbrains.org/groups/biolabs/tools/span-peak-analyzer>"; exit 1;
fi
WORK_DIR=$2
GENOME=$3
CHROM_SIZES=$4
BIN=$5
FDR=$6
if [[ -z "$FDR" ]]; then
    FDR=0.05
fi
GAP=$7
if [[ -z "$GAP" ]]; then
    GAP=5
fi

cd ${WORK_DIR}

TASKS=()
for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    INPUT=$(python ./util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE}: control file: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${FDR}_${GAP}

    if [[ ! -f ${ID}.peak ]]; then
        if [[ -f "${INPUT}" ]]; then
            echo "${FILE}: control file found: ${INPUT}"
            java -Xmx8G -jar ${SPAN_JAR_PATH} analyze -t ${FILE} -c ${INPUT} --chrom.sizes ${CHROM_SIZES} \
                --bin ${BIN} --fdr ${FDR} --gap ${GAP} \
                --peaks ${ID}.peak \
                --threads 4 2&>1 |\ tee ${NAME}_span_${GENOME}.log
        else
            echo "${FILE}: no control file"
            java -jar ${SPAN_JAR_PATH} analyze -t ${FILE} --chrom.sizes ${CHROM_SIZES} \
                --bin ${BIN} --fdr ${FDR} --gap ${GAP} \
                --peaks ${ID}.peak \
                --threads 4 2&>1 |\ tee ${NAME}_span_${GENOME}.log
        fi
    fi
done

check_logs

>&2 echo "Done. Batch SPAN $@"