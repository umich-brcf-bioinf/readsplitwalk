#! /usr/bin/env python

"""
cluster_gaps.py
7/10/2013- cgates/pulintz
7/14/2013 - cgates: refactored to use a Gap object
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

class Gap():

	@staticmethod
	def header():
		return "#chromosome|gap_start|gap_end|gap_width|read_start|read_end|read_width|split_read_name|original_read_name"

	_name_re = re.compile(r"(.+)-([LR])-([\d]+)$")

	def __init__(self, split_read_name, chromosome, read_start, gap_start, gap_end, read_end):
		self._chromosome = chromosome
		self._gap_start = gap_start
		self._gap_end = gap_end
		self._read_start = read_start
		self._read_end = read_end
		self._split_read_name = split_read_name

	def __eq__(self, other):
		if type(other) is type(self):
			return self.__dict__ == other.__dict__
		return False
		
	def _read_width(self):
		return self._read_end - self._read_start
	
	def _gap_width(self):
		return self._gap_end - self._gap_start

	def _original_read_name(self):
		m = self._name_re.match(self._split_read_name)
		(original_read_name, side, split_len) = (m.group(1), m.group(2), m.group(3))
		return original_read_name
	
	def _format(self, delimiter):
		return delimiter.join(
			[self._chromosome, str(self._gap_start), str(self._gap_end), str(self._gap_width()), 
				str(self._read_start), str(self._read_end), str(self._read_width()), self._split_read_name, self._original_read_name()])

class GapUtility():
	
	def __init__(self, original_read_len, delimiter):
		self._original_read_len = int(original_read_len)
		self._delimiter = delimiter

	def build_gap(self, sam_line):
		(split_read_name, flag, transcript_name, start_pos, mapq_score, cigar, rnext, pnext, tlen, seq) = sam_line.split(self._delimiter)[0:10]
		l_gap_pos = int(start_pos) + len(seq)
		r_gap_pos = int(pnext)
		leftmost_start = int(start_pos)
		rightmost_end = int(pnext) + (self._original_read_len - len(seq))
		return Gap(split_read_name, transcript_name, leftmost_start, l_gap_pos, r_gap_pos, rightmost_end)

	def samfile_to_gaps(self, sam_file):
		gaps = []
		for line in sam_file:
			if line.startswith("@"):
				continue
			tlen = int(line.split(self._delimiter)[8])
			if tlen > 0: 
				gaps.append(self.build_gap(line))
		return gaps

	def sort_gaps(self, gaps):
		return sorted(gaps, key=lambda gap: (gap._chromosome, gap._gap_start))

	def write_gap_file(self, sorted_gaps, writer, additional_header_lines=[]):
		for line in additional_header_lines:
			writer.write("#")
			writer.write(line)
			writer.write("\n")
		writer.write(Gap.header())
		writer.write("\n")
		for gap in sorted_gaps:
			writer.write(gap._format(self._delimiter))
			writer.write("\n")

def main(sam_file_name, original_read_len, gap_file_name, delimiter):
 	logger = StdErrLogger(verbose=True)
	logger.log(" ".join(sys.argv), verbose=False)
	header_lines = [str(datetime.datetime.today()), " ".join(sys.argv)] 
	gap_utility = GapUtility(original_read_len, delimiter)
	
	logger.log("parsing sam file")
	sam_file = open(sam_file_name,"r")
	gaps = gap_utility.samfile_to_gaps(sam_file)
	sam_file.close()

	logger.log("sorting {0} gaps".format(len(gaps)))	
	gaps = gap_utility.sort_gaps(gaps)

	logger.log("writing {0} gaps to file".format(len(gaps)))
	gap_file = open(gap_file_name, "w")
	gap_utility.write_gap_file(gaps, gap_file, header_lines)
	gap_file.close()
	
	
	logger.log("complete")
	
if __name__ == "__main__":
	basename = os.path.basename(sys.argv[0])
	if (len(sys.argv) != 4):
		print "usage: {0} [sam_file] [original_read_len] [gap_file]".format(basename)
		sys.exit() 

	(sam_file_name, original_read_len, gap_file_name) = sys.argv[1:]
	sam_file_name = os.path.abspath(sam_file_name)
	gap_file_name = os.path.abspath(gap_file_name)

	main(sam_file_name, original_read_len, gap_file_name, "\t") 
	print basename + " done."
