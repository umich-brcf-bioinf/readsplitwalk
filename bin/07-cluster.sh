#!/bin/bash

BIN_DIR="`dirname $0`"
SCRIPT_NAME="`basename $0 .sh`"
ORIGINAL_READ_LEN=$1
PREDECSSOR_SUFFIX=.06-postprocess.sam
mkdir -p logs
chmod g+rw logs
LOG_FILE=logs/${SCRIPT_NAME}.log

function cluster {
	SAM=$1
	OUTPUT_BASE_NAME=`basename $1 $PREDECSSOR_SUFFIX`.${SCRIPT_NAME}
	echo clustering $SAM ...
	"${BIN_DIR}"/cluster_gaps.py ${SAM} ${ORIGINAL_READ_LEN} ${OUTPUT_BASE_NAME}.tab ${OUTPUT_BASE_NAME}.sam
	
	HEADER="#transcript	cluster	sample_name	gap_start_mean	gap_start_stdev	gap_length_min	gap_length_mean	gap_length_max	gap_length_stdev	split_read_count	original_read_count" 
	echo "${HEADER}" > ${OUTPUT_BASE_NAME}.sample-groups.tab
    #tab file header
	#chromosome      cluster	sample_name    gap_start       gap_end gap_width       read_start      read_end        read_width      split_read_name original_read_name      
    
	echo "sorting ${OUTPUT_BASE_NAME}.tab ..."
	#replace with cat with pv when available
	time cat ${OUTPUT_BASE_NAME}.tab | sort -k1,1 -k2,2n -k3,3 > ${OUTPUT_BASE_NAME}.sorted.tab
	
	echo "building ${OUTPUT_BASE_NAME}.sample-groups.tab ..."
	time bedtools groupby -i ${OUTPUT_BASE_NAME}.sorted.tab \
        -grp 1,2,3 \
        -opCols 4,4,6,6,6,6,11,11 \
        -ops mean,stdev,min,mean,max,stdev,count,count_distinct \
        >> ${OUTPUT_BASE_NAME}.sample-groups.tab

	echo "building ${OUTPUT_BASE_NAME}.cluster-groups.tab ..."
	HEADER="#transcript	cluster	distinct_sample_count	gap_start_mean	gap_start_stdev	gap_length_min	gap_length_mean	gap_length_max	gap_length_stdev	split_read_count	original_read_count" 
	echo "${HEADER}" > ${OUTPUT_BASE_NAME}.cluster-groups.tab
	time bedtools groupby -i ${OUTPUT_BASE_NAME}..sorted.tab \
		-grp 1,2 \
		-opCols 3,4,4,6,6,6,6,11,11 \
		-ops count_distinct,mean,stdev,min,mean,max,stdev,count,count_distinct \
		>> ${OUTPUT_BASE_NAME}.cluster-groups.tab

}

function convert_sam_to_bam {
	BASENAME=`basename $1 .sam`
	BAMNAME=${BASENAME}.bam
	SORTEDNAME=${BASENAME}.sorted
	INDEXEDNAME=${SORTEDNAME}.bam
	echo Converting $1 to ${BAMNAME}...
	time samtools view -bS $1 > ${BAMNAME}
	time samtools sort ${BAMNAME} ${SORTEDNAME}
	time samtools index ${INDEXEDNAME}
}


(
echo teeing to $LOG_FILE
date

cluster all_samples_merged.06-postprocess.sam

convert_sam_to_bam all_samples_merged.${SCRIPT_NAME}.sam

chmod g+rw ${LOG_FILE} *.${SCRIPT_NAME}.* 
date
echo done.
) 2>&1 | tee ${LOG_FILE}

