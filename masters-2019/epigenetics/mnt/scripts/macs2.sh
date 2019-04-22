#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

which macs2 &>/dev/null || { echo "ERROR: MACS2 not found! Download MACS2: <https://github.com/taoliu/MACS/wiki/Install-macs2>"; exit 1; }
which Rscript &>/dev/null || { echo "ERROR: R not found!"; exit 1; }

# Load utils
source $(dirname $0)/util.sh

>&2 echo "Batch macs2 $@"
if [[ $# -lt 4 ]]; then
    echo "Need 4 parameters! <work_dir> <genome> <suffix> <params_str>"
    exit 1
fi

WORK_DIR=$1
GENOME=$2
SUFFIX=$3
PARAMS=$4

SPECIES=$(python $(dirname $0)/util.py macs_species ${GENOME})

cd ${WORK_DIR}

for FILE in $(find . -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
    INPUT=$(python $(dirname $0)/util.py find_input ${WORK_DIR}/${FILE})
    echo "${FILE}: control file: ${INPUT}"

    NAME=${FILE%%.bam} # file name without extension
    ID=${NAME}_${SUFFIX}

    PEAKS_FILE=$(find . -name "${ID}*.*Peak")
    if [[ -z ${PEAKS_FILE} ]]; then
        TMP_DIR=$(mktemp -d ${WORK_DIR}/macs2.XXXXXXXX)
        mkdir -p ${TMP_DIR}

        echo "Macs2 TMP_DIR: ${TMP_DIR}"
        if [[ -f "${INPUT}" ]]; then
            echo "${FILE}: control file found: ${INPUT}"
            macs2 callpeak --tempdir ${TMP_DIR} -t ${FILE} -c ${INPUT} -f BAM -g ${SPECIES} -n ${ID} ${PARAMS} 2>&1 |\
                tee ${ID}_macs2.log
        else
            echo "${FILE}: no control file"
            macs2 callpeak --tempdir ${TMP_DIR} -t ${FILE} -f BAM -g ${SPECIES} -n ${ID} ${PARAMS} 2>&1 |\
                tee ${ID}_macs2.log
        fi
        # Cleanup
        rm -rf ${TMP_DIR}
     fi
done

check_logs

# Create pdf reports
MODELS=$(ls *.r); for M in ${MODELS[@]}; do Rscript $M; done

>&2 echo "Done. Batch macs2 $@"