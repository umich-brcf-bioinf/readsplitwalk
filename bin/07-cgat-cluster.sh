#!/bin/bash

BIN_DIR="`dirname $0`"
SCRIPT_NAME="`basename $0 .sh`"
ORIGINAL_READ_LEN=$1
SUFFIX=06-postprocess.sam

function cluster {
	SAM=$1
	OUTPUT_BASE_NAME=${2}.${SCRIPT_NAME}
	echo Clustering $SAM ...
	"${BIN_DIR}"/cluster_gaps.py ${SAM} ${ORIGINAL_READ_LEN} ${OUTPUT_BASE_NAME}.gaps.tmp
	bedtools cluster -i ${OUTPUT_BASE_NAME}.gaps.tmp > ${OUTPUT_BASE_NAME}.cluster.tmp
	bedtools groupby -i ${OUTPUT_BASE_NAME}.cluster.tmp -grp 10 -opCols 1,2,3,4,4,5,6,7,8,9 -ops distinct,min,max,mean,stdev,min,max,mean,count,count_distinct > ${OUTPUT_BASE_NAME}.groupby.tmp
	awk '{print $0,"\t",(($4-$3)-$5)}' ${OUTPUT_BASE_NAME}.groupby.tmp > ${OUTPUT_BASE_NAME}.tsv
	
#	rm ${OUTPUT_BASE_NAME}.*.tmp
}


cluster Sample_21786_ALL.${SUFFIX} Sample_21786_ALL
cluster Sample_21787_ALL.${SUFFIX} Sample_21787_ALL
exit 1
cluster Sample_21788_ALL.${SUFFIX} Sample_21788_ALL
cluster Sample_21789_ALL.${SUFFIX} Sample_21789_ALL
cluster Sample_21790_ALL.${SUFFIX} Sample_21790_ALL
cluster Sample_21791_ALL.${SUFFIX} Sample_21791_ALL
cluster Sample_21792_ALL.${SUFFIX} Sample_21792_ALL
cluster Sample_21793_ALL.${SUFFIX} Sample_21793_ALL
cluster Sample_21794_ALL.${SUFFIX} Sample_21794_ALL
cluster Sample_21795_ALL.${SUFFIX} Sample_21795_ALL
cluster Sample_21796_ALL.${SUFFIX} Sample_21796_ALL
cluster Sample_21797_ALL.${SUFFIX} Sample_21797_ALL

echo Done.
