#! /usr/bin/env python

"""
cluster_gaps.py
7/10/2013 - cgates/pulintz
"""

import datetime
import os
import re
import resource
import sys
import traceback

class StdErrLogger():
        """Writes basic utilization data to stderr"""
        def __init__(self, verbose = False):
                self._verbose = verbose

        def _log(self, message):
                print >> sys.stderr, "{0}|{1}".format(datetime.datetime.today(), message)


        def log(self, message, verbose=None):
		verbose = self._verbose if verbose is None else verbose 
                if (verbose):
                        usage = resource.getrusage(resource.RUSAGE_SELF)
                        memory_used = usage.ru_maxrss/1024
                        function_name = traceback.extract_stack()[-2:-1][0][2]
                        self._log("usertime(s)={0:.0f}|systime(s)={1:.0f}|peak_memory_used(mb)={2}|{3}|{4}".format(usage.ru_utime, usage.ru_stime, memory_used, function_name, message))
                else:
                        self._log(message)

class BedUtility():
	def __init__(self, original_read_len, delimiter):
		self._original_read_len = int(original_read_len)
		self._delimiter = delimiter
		self._name_re = re.compile(r"(.+)-([LR])-([\d]+)$")

	def bed_line(self, sam_line):
		(split_read_name, flag, transcript_name, start_pos, mapq_score, cigar, rnext, pnext, tlen, seq) = sam_line.split(self._delimiter)[0:10]
		l_gap_pos = int(start_pos) + len(seq)
		r_gap_pos = int(pnext)
		leftmost_start = int(start_pos)
		rightmost_end = int(pnext) + (self._original_read_len - len(seq))
		strand = "+" if tlen > 0 else "-"
		m = self._name_re.match(split_read_name)
                (original_read_name, side, split_len) = (m.group(1), m.group(2), m.group(3))
		return self._delimiter.join([transcript_name, str(l_gap_pos), str(r_gap_pos), split_read_name, mapq_score, strand, str(leftmost_start), str(rightmost_end), original_read_name])

	def samfile_to_bed(self, sam_file):
		bed_lines=["#chrom|gapStart|gapEnd|readName|score|strand|readStart|readEnd"]
		for line in sam_file:
			if line.startswith("@"):
				continue
			tlen = int(line.split(self._delimiter)[8])
			if tlen > 0: 
				bed_lines.append(self.bed_line(line))
		return bed_lines

	def sort_bed_lines(self, bed_lines):
		def line_key(line):
			if line.startswith("#"):
				return ("", 0) 
			(chrom, start) = line.split(self._delimiter)[0:2]
			return (chrom, int(start))
		return sorted(bed_lines, key=line_key)

	def write_bed_file(self, bed_lines, writer, additional_header_lines=[]):
		for line in additional_header_lines:
			writer.write("#")
			writer.write(line)
			writer.write("\n")
		for line in bed_lines:
			writer.write(line)
			writer.write("\n")

def main(sam_file_name, original_read_len, bed_file_name, delimiter):
 	logger = StdErrLogger(verbose=True)
	logger.log(" ".join(sys.argv), verbose=False)
	header_lines = [str(datetime.datetime.today()), " ".join(sys.argv)] 
	bed_utility = BedUtility(original_read_len, delimiter)
	
	logger.log("parsing sam file")
	sam_file = open(sam_file_name,"r")
	bed_lines = bed_utility.samfile_to_bed(sam_file)
	sam_file.close()

	logger.log("sorting bed lines")	
	bed_lines = bed_utility.sort_bed_lines(bed_lines)

	logger.log("writing bed file")
	bed_file = open(bed_file_name, "w")
	bed_utility.write_bed_file(bed_lines, bed_file, header_lines)
	bed_file.close()
	
	logger.log("complete")
	
if __name__ == "__main__":
	if (len(sys.argv) != 4):
		print "usage: {0} [sam_file] [original_read_len] [bed_file]".format(os.path.basename(sys.argv[0]))
		sys.exit() 

	(sam_file_name, original_read_len, bed_file_name) = sys.argv[1:]
	sam_file_name = os.path.abspath(sam_file_name)
	bed_file_name = os.path.abspath(bed_file_name)

	main(sam_file_name, original_read_len, bed_file_name, "\t") 
	print "done."
