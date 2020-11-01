#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author os@jetbrains.com.com

# Stop exec on error.
set -e

which java &>/dev/null || { echo "ERROR: java not found!"; exit 1; }

>&2 echo "Batch SPAN $@"
if [[ $# -lt 4 ]]; then
    echo "Need at least 4 parameters! <SPAN_JAR_PATH> <WORK_DIR> <GENOME> <CHROM_SIZES> [<BIN>] [<FDR>] [<GAP>]"
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
if [[ -z "$BIN" ]]; then
    BIN=200
fi
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
    INPUT=$(python $(dirname $0)/util.py find_input ${WORK_DIR}/${FILE})
    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${FDR}_${GAP}

    if [[ ! -f ${ID}.peak ]]; then
        if [[ -f "${INPUT}" ]]; then
            echo "${FILE}: control file found: ${INPUT}"
            java -Xmx16G -jar ${SPAN_JAR_PATH} analyze -t ${FILE} -c ${INPUT} --chrom.sizes ${CHROM_SIZES} \
                --bin ${BIN} --fdr ${FDR} --gap ${GAP} \
                --peaks ${ID}.peak \
                --threads 6 2>&1 | tee ${NAME}_span_${GENOME}.log
        else
            echo "${FILE}: no control file"
            java -Xmx16G -jar ${SPAN_JAR_PATH} analyze -t ${FILE} --chrom.sizes ${CHROM_SIZES} \
                --bin ${BIN} --fdr ${FDR} --gap ${GAP} \
                --peaks ${ID}.peak \
                --threads 6 2>&1 | tee ${NAME}_span_${GENOME}.log
        fi
    fi
done

>&2 echo "Done. Batch SPAN $@"