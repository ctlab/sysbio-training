#!/usr/bin/env bash

WORK_DIR=$(pwd)
GENOME="hg19"
INDEX_FOLDER=/mnt/${GENOME}

echo "Get genome information"
bash ./get_genome.sh ${GENOME} ${INDEX_FOLDER}

echo "Build bowtie indexes"
bash ./index_bowtie.sh ${GENOME} ${INDEX_FOLDER}

echo "QC for reads file"
bash ./fastqc.sh ${WORK_DIR}

echo "Align"
bash ./bowtie.sh ${GENOME} ${INDEX_FOLDER} 0 ${WORK_DIR}

echo "QC"
bash ./fastqc.sh ${WORK_DIR}
bash ./bam_qc.sh /opt/phantompeakqualtools ${WORK_DIR}

echo "Visualization"
bash ./bigwig.sh ${WORK_DIR}

echo "Remove duplicates"
bash ./remove_duplicates.sh /opt/picard.jar ${WORK_DIR}

echo "Peak calling"
bash ./macs2.sh ${WORK_DIR} ${GENOME} "q0.05" "-q 0.05"
bash ./macs2.sh ${WORK_DIR} ${GENOME} "broad_0.1" "--broad --broad-cutoff 0.1"
bash ./sicer.sh ${WORK_DIR} ${GENOME} ${INDEX_FOLDER}/${GENOME}.chrom.sizes 0.05
bash ./span.sh /opt/span.jar ${WORK_DIR} ${GENOME} ${INDEX_FOLDER}/${GENOME}.chrom.sizes 200 0.05 5

echo "DONE"

