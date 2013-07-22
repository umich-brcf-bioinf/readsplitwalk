#!/bin/bash

function merge_sam {
	READ_1=$1
	READ_2=$2
	MERGE_FILE=$3
	echo Merging to ${MERGE_FILE}...
	cp ${READ_1} ${MERGE_FILE}
	grep '^[^@]' ${READ_2} >> ${MERGE_FILE}
}

function samtools_process {
	BASENAME=`basename $1 .sam`
	BAMNAME=${BASENAME}.bam
	SORTEDNAME=${BASENAME}.sorted
	INDEXEDNAME=${SORTEDNAME}.bam
	echo Converting $1 to ${BAMNAME}...
	samtools view -bS $1 > ${BAMNAME}
	samtools sort ${BAMNAME} ${SORTEDNAME}
	samtools index ${INDEXEDNAME}
}

#merge_sam Sample_21786_R1.05-identify_pairs_transcriptome.sam Sample_21786_R2.05-identify_pairs_transcriptome.sam Sample_21786_ALL.06-postprocess.sam
merge_sam Sample_21787_R1.05-identify_pairs_transcriptome.sam Sample_21787_R2.05-identify_pairs_transcriptome.sam Sample_21787_ALL.06-postprocess.sam
#merge_sam Sample_21788_R1.05-identify_pairs_transcriptome.sam Sample_21788_R2.05-identify_pairs_transcriptome.sam Sample_21788_ALL.06-postprocess.sam
#merge_sam Sample_21789_R1.05-identify_pairs_transcriptome.sam Sample_21789_R2.05-identify_pairs_transcriptome.sam Sample_21789_ALL.06-postprocess.sam
#merge_sam Sample_21790_R1.05-identify_pairs_transcriptome.sam Sample_21790_R2.05-identify_pairs_transcriptome.sam Sample_21790_ALL.06-postprocess.sam
# merge_sam Sample_21791_R1.05-identify_pairs_transcriptome.sam Sample_21791_R2.05-identify_pairs_transcriptome.sam Sample_21791_ALL.06-postprocess.sam
# merge_sam Sample_21792_R1.05-identify_pairs_transcriptome.sam Sample_21792_R2.05-identify_pairs_transcriptome.sam Sample_21792_ALL.06-postprocess.sam
# merge_sam Sample_21793_R1.05-identify_pairs_transcriptome.sam Sample_21793_R2.05-identify_pairs_transcriptome.sam Sample_21793_ALL.06-postprocess.sam
# merge_sam Sample_21794_R1.05-identify_pairs_transcriptome.sam Sample_21794_R2.05-identify_pairs_transcriptome.sam Sample_21794_ALL.06-postprocess.sam
# merge_sam Sample_21795_R1.05-identify_pairs_transcriptome.sam Sample_21795_R2.05-identify_pairs_transcriptome.sam Sample_21795_ALL.06-postprocess.sam
# merge_sam Sample_21796_R1.05-identify_pairs_transcriptome.sam Sample_21796_R2.05-identify_pairs_transcriptome.sam Sample_21796_ALL.06-postprocess.sam
merge_sam Sample_21797_R1.05-identify_pairs_transcriptome.sam Sample_21797_R2.05-identify_pairs_transcriptome.sam Sample_21797_ALL.06-postprocess.sam

#samtools_process Sample_21786_ALL.06-postprocess.sam
samtools_process Sample_21787_ALL.06-postprocess.sam
# samtools_process Sample_21788_ALL.06-postprocess.sam
# samtools_process Sample_21789_ALL.06-postprocess.sam
# samtools_process Sample_21790_ALL.06-postprocess.sam
# samtools_process Sample_21791_ALL.06-postprocess.sam
# samtools_process Sample_21792_ALL.06-postprocess.sam
# samtools_process Sample_21793_ALL.06-postprocess.sam
# samtools_process Sample_21794_ALL.06-postprocess.sam
# samtools_process Sample_21795_ALL.06-postprocess.sam
# samtools_process Sample_21796_ALL.06-postprocess.sam
samtools_process Sample_21797_ALL.06-postprocess.sam

echo Done.
