#! /usr/bin/env python

"""
cluster_gaps.py
7/10/2013 - cgates/pulintz
"""

import sys
import os

def bed_line(sam_line, original_read_len, delimiter="\t", name_delimiter="|"):
	(rname, flag, transcript_name, start_pos, mapq, cigar, rnext, pnext, tlen, seq) = sam_line.split(delimiter)[0:10]
	l_gap_pos = int(start_pos) + len(seq)
	r_gap_pos = int(pnext) - 1
	leftmost_start = int(start_pos)
	rightmost_end = int(pnext) + (original_read_len - len(seq))
	name = name_delimiter.join([rname, str(leftmost_start), str(rightmost_end)])  
	return delimiter.join([transcript_name, str(l_gap_pos), str(r_gap_pos), name]) + "\n"


def samfile_to_bedfile(sam_file, original_read_len, writer, field_delimiter="\t", name_delimiter="|"):
	for line in sam_file:
		if line.startswith("@"):
			continue
		tlen = int(line.split(field_delimiter)[8])
		if tlen > 0: 
			writer.write(bed_line(line, original_read_len, field_delimiter, name_delimiter))

def main(sam_file_name, bed_file_name, delimiter="\t"):
	sam_file = open(sam_file_name,"r")
	bed_file = open(bed_file_name, "w")
	
	samfile_to_bedfile(sam_file, bed_file, delimiter)

	sam_file.close()
	bed_file.close()
	
if __name__ == "__main__":
	if (len(sys.argv) != 3):
		print "usage: {0} [sam_file] [bed_file]".format(os.path.basename(sys.argv[0]))
		sys.exit() 

	(sam_file_name, bed_file_name) = sys.argv[1:]
	sam_file_name = os.path.abspath(sam_file_name)
	bed_file_name = os.path.abspath(bed_file_name)

	main(sam_file_name, bed_file_name, "\t", "|") 
	print "done."
