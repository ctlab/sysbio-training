#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

# Checks for errors in logs, stops the world
check_logs()
{
    # IGNORE MACS2 ValueError
    # See for details: https://github.com/JetBrains-Research/washu/issues/14
    # Also ignore SPP failure
    ERRORS=$(find . -name "*.log" |\
        xargs grep -i -E "error|exception|No such file or directory" |\
        grep -v -E "ValueError|WARNING")

    if [[ ! -z "$ERRORS" ]]; then
        echo "ERRORS found"
        echo "$ERRORS"
        exit 1
    fi
}

# Computes and returns pileup file for given BAM file
function pileup(){
    if [[ ! $# -eq 1 ]]; then
        echo "Need 1 argument! <bam_file>"
        exit 1
    fi
    BAM=$1
    PILEUP_DIR="$(dirname ${BAM})/pileup"
    if [[ ! -d ${PILEUP_DIR} ]]; then
        >&2 echo "Create pileup dir ${PILEUP_DIR}"
        mkdir -p ${PILEUP_DIR}
    fi
    NAME=$(basename ${BAM/.bam/_pileup.bed})
    RESULT=${PILEUP_DIR}/${NAME}
    if [[ ! -f ${RESULT} ]]; then
        PILEUP_TMP=$(mktemp pileup.XXXXXX.bed)
        >&2 echo "Calculate ${BAM} pileup file in tmp file: ${PILEUP_TMP}"
        bedtools bamtobed -i ${BAM} > ${PILEUP_TMP}
        # Check that we are the first in async calls, not 100% safe
        if [[ ! -f ${RESULT} ]]; then
            mv ${PILEUP_TMP} ${RESULT}
        else
            >&2 echo "Ignore result, file has been already calculated: ${RESULT}"
        fi
    fi
    echo "${RESULT}"
}