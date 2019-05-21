#!/usr/bin/env bash

# Stop exec on error.
set -e

>&2 echo "pipeline $@"
if [[ $# -lt 5 ]]; then
    echo "Need 5 parameters! <WORK_DIR> <INDEX_FOLDER> <GENOME> <SPAN_JAR> <PICARD_TOOLS_JAR>"
    exit 1
fi

SCRIPTS_DIR=$(dirname $0)

WORK_DIR=$1
INDEX_FOLDER=$2
GENOME=$3
SPAN_JAR=$4
PICARD_TOOLS_JAR=$5

GENOME_INDEX=${INDEX_FOLDER}/${GENOME}
if [[ ! -d ${GENOME_INDEX} ]]; then
    mkdir -p ${GENOME_INDEX}
fi

#echo "Get genome information"
#bash ${SCRIPTS_DIR}/get_genome.sh ${GENOME} ${GENOME_INDEX}
#echo "Build bowtie indexes"
#bash ${SCRIPTS_DIR}/index_bowtie.sh ${GENOME} ${GENOME_INDEX}

echo "QC for reads file"
bash ${SCRIPTS_DIR}/fastqc.sh ${WORK_DIR}

echo "Processing multiqc"
multiqc -f -o ${WORK_DIR}/fastqc ${WORK_DIR}/fastqc

echo "Align"
bash ${SCRIPTS_DIR}/bowtie.sh ${GENOME} ${GENOME_INDEX} 0 ${WORK_DIR}

echo "Processing multiqc"
multiqc -f -o ${WORK_DIR}/bams_qc ${WORK_DIR}/*_bowtie*.log

echo "Visualization"
bash ${SCRIPTS_DIR}/bigwig.sh ${WORK_DIR}

echo "Peak calling"
bash ${SCRIPTS_DIR}/macs2.sh ${WORK_DIR} ${GENOME} "q0.05" "-q 0.05"
bash ${SCRIPTS_DIR}/macs2.sh ${WORK_DIR} ${GENOME} "broad_0.1" "--broad --broad-cutoff 0.1"
bash ${SCRIPTS_DIR}/sicer.sh ${WORK_DIR} ${GENOME} ${GENOME_INDEX}/${GENOME}.chrom.sizes 0.05
bash ${SCRIPTS_DIR}/span.sh ${SPAN_JAR} ${WORK_DIR} ${GENOME} ${GENOME_INDEX}/${GENOME}.chrom.sizes 200 0.05 5

echo "Removing duplicates"
bash ${SCRIPTS_DIR}/remove_duplicates.sh ${PICARD_TOOLS_JAR} ${WORK_DIR}

echo "DONE"

