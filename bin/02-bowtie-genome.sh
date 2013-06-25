#!/bin/bash

SAMPLE_NAME=$1


REFERENCE_INDEX_LOCATION=/ccmb/CoreBA/BioinfCore/Common/DATA/BowtieData/mm9/bowtie1/mm9
MAX_DISTANCE=39999
RSW_HOME=`dirname $0`/../
SCRIPT_NAME=`basename $0 .sh`
mkdir -p logs

INPUT_FILE=${SAMPLE_NAME}.01-bowtie-transcriptome.fastq
OUTPUT_FILE=${SAMPLE_NAME}.${SCRIPT_NAME}
LOG_FILE=logs/${SAMPLE_NAME}.${SCRIPT_NAME}.log

if [ ! -f ${INPUT_FILE} ]; then
	echo Input file [${INPUT_FILE}] not found. Check specified sample name [${SAMPLE_NAME}] is correct.
	exit 1
fi

(
echo $0 $SAMPLE_NAME
echo teeing to $LOG_FILE
date
time bowtie -t -v 0 -k 11 -m 10 --best ${REFERENCE_INDEX_LOCATION} -q ${INPUT_FILE} --un ${OUTPUT_FILE}.fastq ${OUTPUT_FILE}.tmp.bowtie
chmod g+rw ${LOG_FILE} ${OUTPUT_FILE}*
date
echo done.
) 2>&1 | tee ${LOG_FILE}
