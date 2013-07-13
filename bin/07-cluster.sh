#!/bin/bash

ORIGINAL_READ_LEN=$1

function cluster {
	SAM=$1
	OUTPUT_BASE_NAME=$2
	echo Clustering $SAM ...
	/ccmb/CoreBA/BioinfCore/Common/DATA/Justin_Kaufman_Data/readSplitWalk/bin/cluster_gaps.py ${SAM} ${ORIGINAL_READ_LEN} ${OUTPUT_BASE_NAME}.gaps.tmp
	bedtools cluster -i ${OUTPUT_BASE_NAME}.gaps.tmp > ${OUTPUT_BASE_NAME}.cluster.tmp
	bedtools groupby -i ${OUTPUT_BASE_NAME}.cluster.tmp -grp 10 -opCols 1,2,3,4,4,5,6,7,8,9 -ops distinct,min,max,mean,stdev,min,max,mean,count,count_distinct > ${OUTPUT_BASE_NAME}.groupby.tmp
	awk '{print $0,"\t",(($4-$3)-$5)}' ${OUTPUT_BASE_NAME}.groupby.tmp > ${OUTPUT_BASE_NAME}.tsv
	
	rm ${OUTPUT_BASE_NAME}.*.tmp
}


cluster Sample_21786_ALL.06-postprocess.sam Sample_21786_ALL.07-cluster
cluster Sample_21787_ALL.06-postprocess.sam Sample_21787_ALL.07-cluster
cluster Sample_21788_ALL.06-postprocess.sam Sample_21788_ALL.07-cluster
cluster Sample_21789_ALL.06-postprocess.sam Sample_21789_ALL.07-cluster
cluster Sample_21790_ALL.06-postprocess.sam Sample_21790_ALL.07-cluster
cluster Sample_21791_ALL.06-postprocess.sam Sample_21791_ALL.07-cluster
cluster Sample_21792_ALL.06-postprocess.sam Sample_21792_ALL.07-cluster
cluster Sample_21793_ALL.06-postprocess.sam Sample_21793_ALL.07-cluster
cluster Sample_21794_ALL.06-postprocess.sam Sample_21794_ALL.07-cluster
cluster Sample_21795_ALL.06-postprocess.sam Sample_21795_ALL.07-cluster
cluster Sample_21796_ALL.06-postprocess.sam Sample_21796_ALL.07-cluster
cluster Sample_21797_ALL.06-postprocess.sam Sample_21797_ALL.07-cluster

echo Done.
