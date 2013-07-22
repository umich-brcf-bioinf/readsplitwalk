#!/bin/bash

SAMPLE_NAME=$1
SPLIT_MARGIN=$2

RSW_HOME=`dirname $0`/../
SCRIPT_NAME=`basename $0 .sh`
mkdir -p logs

INPUT_FILE=${SAMPLE_NAME}.02-bowtie-genome.fastq
OUTPUT_FILE=${SAMPLE_NAME}.${SCRIPT_NAME}.fastq
LOG_FILE=logs/${SAMPLE_NAME}.${SCRIPT_NAME}.log

if [ ! -f ${INPUT_FILE} ]; then
	echo Input file [${INPUT_FILE}] not found. Check specified sample name [${SAMPLE_NAME}] is correct.
	exit 1
fi

(
echo $0 $SAMPLE_NAME
echo teeing to $LOG_FILE
date
time python ${RSW_HOME}/bin/split_read.py ${INPUT_FILE} ${OUTPUT_FILE} ${SPLIT_MARGIN}
chmod g+rw ${LOG_FILE} ${OUTPUT_FILE}
date
echo done.
) 2>&1 | tee ${LOG_FILE}
