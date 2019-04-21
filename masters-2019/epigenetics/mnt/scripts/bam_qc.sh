#!/usr/bin/env bash
# Original https://github.com/JetBrains-Research/washu
# Author oleg.shpynov@jetbrains.com

which bedtools &>/dev/null || { echo "ERROR: bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which samtools &>/dev/null || { echo "ERROR: samtools not found! Download samtools: <http://www.htslib.org/doc/samtools.html>"; exit 1; }
which Rscript &>/dev/null || { echo "ERROR: R not found!"; exit 1; }

>&2 echo "bam_qc: $@"

# Load utils
source ./util.sh

>&2 echo "Batch bam_qc $@"
if [[ $# -lt 2 ]]; then
    echo "Need 2 parameters! <phantompeakqualtools> <work_dir>"
    exit 1
fi

PHANTOMPEAKQUALTOOLS=$1
WORK_DIRS=$2

for FILE in $(find . -name '*.bam' | sed 's#\./##g')
do :
    NAME=${FILE%%.bam} # file name without extension

    export TMP_DIR=$(type job_tmp_dir &>/dev/null && echo "\$(tmp_dir)" || echo "/tmp")

    #col.	abbreviation	description
    #1	Filename	tagAlign/BAM filename
    #2	numReads	effective sequencing depth i.e. total number of mapped reads in input file
    #3	estFragLen	comma separated strand cross-correlation peak(s) in decreasing order of correlation.
    #4	corr_estFragLen	comma separated strand cross-correlation value(s) in decreasing order (COL2 follows the same order)
    #5	phantomPeak	Read length/phantom peak strand shift
    #6	corr_phantomPeak	Correlation value at phantom peak
    #7	argmin_corr	strand shift at which cross-correlation is lowest
    #8	min_corr	minimum value of cross-correlation
    #9	NSC	Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8
    #10	RSC	Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8)
    #11	QualityTag	Quality tag based on thresholded RSC (codes= -2:veryLow, -1:Low, 0:Medium, 1:High, 2:veryHigh)
    Rscript ${PHANTOMPEAKQUALTOOLS}/run_spp.R -c=${FILE} -savp -out=${NAME}.phantom.tsv 2&>1 | tee ${NAME}_bam_qc.log
done

check_logs

>&2 echo "Done. Batch bams2_qc $@"