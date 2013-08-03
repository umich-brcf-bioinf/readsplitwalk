#!/bin/bash

BIN_DIR="`dirname $0`"
SCRIPT_NAME="`basename $0 .sh`"
PREDECESSOR_PREFIX=all_samples_merged.07-cluster
OUTPUT_BASE_NAME=all_samples_merged.${SCRIPT_NAME}
TRANSCRIPT_MAPPING_FILE=${BIN_DIR}/../reference_data/ensembl_GRCm38-p1_cDNA_gene_mapping.txt

mkdir -p logs
chmod g+rw logs
LOG_FILE=logs/${SCRIPT_NAME}.log


function add_gene_symbol {
	INPUT_FILE=$1
	OUTPUT_FILE=$2
	TRANSCRIPT_COLUMN_INDEX=0
	echo "Adding gene symbol to ${INPUT_FILE}..."
	
	"${BIN_DIR}"/transcript_to_gene_symbol.py ${INPUT_FILE} ${TRANSCRIPT_COLUMN_INDEX} ${TRANSCRIPT_MAPPING_FILE} ${OUTPUT_FILE}

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

add_gene_symbol ${PREDECESSOR_PREFIX}.sorted.tab ${OUTPUT_BASE_NAME}.tab
add_gene_symbol ${PREDECESSOR_PREFIX}.cluster-groups.tab ${OUTPUT_BASE_NAME}.cluster-groups.tab
add_gene_symbol ${PREDECESSOR_PREFIX}.sample-groups.tab ${OUTPUT_BASE_NAME}.sample-groups.tab

chmod g+rw ${LOG_FILE} *.${SCRIPT_NAME}.* 
date
echo done.
) 2>&1 | tee ${LOG_FILE}

