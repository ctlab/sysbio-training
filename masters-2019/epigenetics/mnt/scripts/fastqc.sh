#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com
# Author roman.chernyatchik@jetbrains.com

which fastqc &>/dev/null || { echo "ERROR: fastqc not found!"; exit 1; }

# Load utils
source $(dirname $0)/util.sh

>&2 echo "Batch fastqc $@"
if [[ $# -lt 1 ]]; then
    echo "Need one parameter! <WORK_DIR>"
    exit 1
fi
WORK_DIR=$1
cd ${WORK_DIR}

RESULTS_DIR="fastqc"
if [[ -d "${RESULTS_DIR}" ]]; then
    echo "[Skipped] ${RESULTS_DIR} was already processed"
    exit 0
else
    mkdir -p "${RESULTS_DIR}"
fi

for FILE in $(find . -name '*.f*q' | sed 's#\./##g')
do :
    FILE_NAME=${FILE##*/}
    NAME=${FILE_NAME%%.f*q} # file name without extension

    fastqc --outdir "${RESULTS_DIR}" "${FILE}" 2>&1 | tee ${NAME}_fastqc.log
done

check_logs

>&2 echo "Done. Batch fastqc $@"
