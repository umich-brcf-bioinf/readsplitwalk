#!/bin/bash

SCRIPTDIR=`dirname $0`

function merge_sam {
	READ_1=$1
	READ_2=$2
	MERGE_FILE=$3
	echo Merging to ${MERGE_FILE}...
	cp ${READ_1} ${MERGE_FILE}
	grep '^[^@]' ${READ_2} >> ${MERGE_FILE}
}

function convert_sam_to_bam {
	BASENAME=`basename $1 .sam`
	BAMNAME=${BASENAME}.bam
	SORTEDNAME=${BASENAME}.sorted
	INDEXEDNAME=${SORTEDNAME}.bam
	echo Converting $1 to ${BAMNAME}...
	samtools view -bS $1 > ${BAMNAME}
	samtools sort ${BAMNAME} ${SORTEDNAME}
	samtools index ${INDEXEDNAME}
}

function add_read_group {
	BASENAME=`basename $1 .sam`
	RGBASENAME=`basename $1 _ALL.06-postprocess.sam`
	RGBASENAMEOUT=${BASENAME}.rg.sam
	echo Adding readgroup info to $1...
	python ${SCRIPTDIR}/add_readgroup_to_sam.py $1 ${RGBASENAME} ${RGBASENAME} ${RGBASENAMEOUT}
}



merge_sam Sample_21786_R1.05-identify_pairs_transcriptome.sam Sample_21786_R2.05-identify_pairs_transcriptome.sam Sample_21786_ALL.06-postprocess.sam
merge_sam Sample_21787_R1.05-identify_pairs_transcriptome.sam Sample_21787_R2.05-identify_pairs_transcriptome.sam Sample_21787_ALL.06-postprocess.sam
merge_sam Sample_21788_R1.05-identify_pairs_transcriptome.sam Sample_21788_R2.05-identify_pairs_transcriptome.sam Sample_21788_ALL.06-postprocess.sam
merge_sam Sample_21789_R1.05-identify_pairs_transcriptome.sam Sample_21789_R2.05-identify_pairs_transcriptome.sam Sample_21789_ALL.06-postprocess.sam
merge_sam Sample_21790_R1.05-identify_pairs_transcriptome.sam Sample_21790_R2.05-identify_pairs_transcriptome.sam Sample_21790_ALL.06-postprocess.sam
merge_sam Sample_21791_R1.05-identify_pairs_transcriptome.sam Sample_21791_R2.05-identify_pairs_transcriptome.sam Sample_21791_ALL.06-postprocess.sam
merge_sam Sample_21792_R1.05-identify_pairs_transcriptome.sam Sample_21792_R2.05-identify_pairs_transcriptome.sam Sample_21792_ALL.06-postprocess.sam
merge_sam Sample_21793_R1.05-identify_pairs_transcriptome.sam Sample_21793_R2.05-identify_pairs_transcriptome.sam Sample_21793_ALL.06-postprocess.sam
merge_sam Sample_21794_R1.05-identify_pairs_transcriptome.sam Sample_21794_R2.05-identify_pairs_transcriptome.sam Sample_21794_ALL.06-postprocess.sam
merge_sam Sample_21795_R1.05-identify_pairs_transcriptome.sam Sample_21795_R2.05-identify_pairs_transcriptome.sam Sample_21795_ALL.06-postprocess.sam
merge_sam Sample_21796_R1.05-identify_pairs_transcriptome.sam Sample_21796_R2.05-identify_pairs_transcriptome.sam Sample_21796_ALL.06-postprocess.sam
merge_sam Sample_21797_R1.05-identify_pairs_transcriptome.sam Sample_21797_R2.05-identify_pairs_transcriptome.sam Sample_21797_ALL.06-postprocess.sam

add_read_group Sample_21786_ALL.06-postprocess.sam
add_read_group Sample_21787_ALL.06-postprocess.sam
add_read_group Sample_21788_ALL.06-postprocess.sam
add_read_group Sample_21789_ALL.06-postprocess.sam
add_read_group Sample_21790_ALL.06-postprocess.sam
add_read_group Sample_21791_ALL.06-postprocess.sam
add_read_group Sample_21792_ALL.06-postprocess.sam
add_read_group Sample_21793_ALL.06-postprocess.sam
add_read_group Sample_21794_ALL.06-postprocess.sam
add_read_group Sample_21795_ALL.06-postprocess.sam
add_read_group Sample_21796_ALL.06-postprocess.sam
add_read_group Sample_21797_ALL.06-postprocess.sam

echo 'Merging all headers...'
grep --no-filename '^@' Sample_*_ALL.06-postprocess.rg.sam | grep -v '^@PG' | sort -k1,1r | uniq > all_samples_merged.06-postprocess.sam
echo 'Merging all alignments...'
grep --no-filename -v '^@' Sample_*_ALL.06-postprocess.rg.sam >> all_samples_merged.06-postprocess.sam

#echo Converting to bam files...
#convert_sam_to_bam all_samples_merged.06-postprocess.sam

echo Done.
