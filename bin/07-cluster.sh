#!/bin/bash

BIN_DIR="`dirname $0`"
SCRIPT_NAME="`basename $0 .sh`"
ORIGINAL_READ_LEN=$1
SUFFIX=06-postprocess.sam
mkdir -p logs
LOG_FILE=logs/${SCRIPT_NAME}.log

function cluster {
	SAM=$1
	SAMPLE_NAME=${2}
	OUTPUT_BASE_NAME=${SAMPLE_NAME}.${SCRIPT_NAME}
	echo ================================
	echo gapping $SAM ...
	"${BIN_DIR}"/cluster_gaps.py ${SAM} ${ORIGINAL_READ_LEN} ${OUTPUT_BASE_NAME}.gaps.tmp
	echo clustering $SAM ...
	#chromosome|gap_start|gap_end|gap_width|read_start|read_end|read_width|split_read_name|original_read_name
	bedtools cluster -i ${OUTPUT_BASE_NAME}.gaps.tmp > ${OUTPUT_BASE_NAME}.cluster.tmp

	echo groupby $SAM ...	
	HEADER="cluster_id	transcript	read_start_min	gap_start_min	gap_start_mean	gap_start_stdev	gap_end_mean	gap_end_stdev	gap_end_max	read_end_max	gap_length_mean	gap_length_stdev	read_length_mean	split_read_count	original_read_count" 
	echo "#${HEADER}" > ${OUTPUT_BASE_NAME}.groupby.tmp
	bedtools groupby -i ${OUTPUT_BASE_NAME}.cluster.tmp -grp 10 -opCols 1,5,2,2,2,3,3,3,6,4,4,7,8,9 -ops distinct,min,min,mean,stdev,mean,stdev,max,max,mean,stdev,mean,count,count_distinct > ${OUTPUT_BASE_NAME}.groupby.tmp
	echo finalizing $SAM ...
	echo "#sample	${HEADER}" > ${OUTPUT_BASE_NAME}.tsv
	awk -v sample=${SAMPLE_NAME} '{print sample,"\t",$0}' ${OUTPUT_BASE_NAME}.groupby.tmp >> ${OUTPUT_BASE_NAME}.tsv
	
	rm ${OUTPUT_BASE_NAME}.*.tmp
}



(
echo teeing to $LOG_FILE
date

cluster Sample_21786_ALL.${SUFFIX} Sample_21786_ALL &
cluster Sample_21787_ALL.${SUFFIX} Sample_21787_ALL &
cluster Sample_21788_ALL.${SUFFIX} Sample_21788_ALL &
cluster Sample_21789_ALL.${SUFFIX} Sample_21789_ALL &
cluster Sample_21790_ALL.${SUFFIX} Sample_21790_ALL &
cluster Sample_21791_ALL.${SUFFIX} Sample_21791_ALL &
cluster Sample_21792_ALL.${SUFFIX} Sample_21792_ALL &
cluster Sample_21793_ALL.${SUFFIX} Sample_21793_ALL &
cluster Sample_21794_ALL.${SUFFIX} Sample_21794_ALL &
cluster Sample_21795_ALL.${SUFFIX} Sample_21795_ALL &
cluster Sample_21796_ALL.${SUFFIX} Sample_21796_ALL &
cluster Sample_21797_ALL.${SUFFIX} Sample_21797_ALL &

wait

awk 'FNR < 2' Sample_21786_ALL.*.tsv > ALL.tsv
awk 'FNR > 1' Sample*.tsv >> ALL.tsv
chmod g+rw ${LOG_FILE} *.${SCRIPT_NAME}.* ALL.tsv
date
echo done.
) 2>&1 | tee ${LOG_FILE}

