#!/bin/bash

SAMPLE_NAME=$1
READ_LENGTH=$2

MIN_DISTANCE=2
MAX_DISTANCE=39999
RSW_HOME=`dirname $0`/../
SCRIPT_NAME=`basename $0 .sh`
mkdir -p logs

INPUT_FILE=${SAMPLE_NAME}.04-bowtie_align_splits_to_transcriptome_1mm.sam
OUTPUT_FILE=${SAMPLE_NAME}.${SCRIPT_NAME}.rsw
LOG_FILE=logs/${SAMPLE_NAME}.${SCRIPT_NAME}.log

if [ ! -f ${INPUT_FILE} ]; then
	echo Input file [${INPUT_FILE}] not found. Check specified sample name [${SAMPLE_NAME}] is correct.
	exit 1
fi

(
echo $0 $SAMPLE_NAME
echo teeing to $LOG_FILE
date
time python ${RSW_HOME}/bin/identify_pairs.py ${INPUT_FILE} ${OUTPUT_FILE} ${READ_LENGTH} ${MIN_DISTANCE} ${MAX_DISTANCE} 
chmod g+rw ${LOG_FILE} ${OUTPUT_FILE}
date
echo done.
) 2>&1 | tee ${LOG_FILE}

