#!/bin/bash

## batch_fastqc.sh <input_dir>
# Runs fastqc program on a directory of fastq files 
# (eg, raw_fastq/filtered_fastq) 
# Must create a fastqc/ directory for the output files
## Ana / October 19, 2011

ARGS=2;
if [ $# -ne "$ARGS" ]
then
  echo "Usage: `basename $0` input_dir zip";
  exit $E_BADARGS;
fi

if [ $2 == "TRUE" ] 
then
  fastq=$1*.fastq.gz;
  fastqc -o fastqc --noextract -f fastq $fastq;
else
  fastq=$1*.fastq;
  fastqc -o fastqc -f fastq $fastq;
fi
