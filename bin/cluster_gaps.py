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

# pylint: disable=R0903
class StdErrLogger():
    """Writes basic utilization data to stderr"""
    def __init__(self, verbose = False):
        self._verbose = verbose

    def log(self, message, verbose=None):
        """Logs message to std err with optional stats on peak utilization"""
        verbose = self._verbose if verbose is None else verbose
        if (verbose):
            usage = resource.getrusage(resource.RUSAGE_SELF)
            memory_used = usage.ru_maxrss/1024
            function_name = traceback.extract_stack()[-2:-1][0][2]
            message = "usertime(s)={0:.0f}|systime(s)={1:.0f}|"\
                "peak_memory_used(mb)={2}|{3}|{4}". \
                format(usage.ru_utime, usage.ru_stime, memory_used, 
                    function_name, message)
        sys.stderr.write("{0}|{1}\n".format(datetime.datetime.today(), message))

class Gap():
    """Models an aligned pair of split reads"""
    
    @staticmethod
    def header():
        # pylint: disable=line-too-long
        return "#chromosome|gap_start|gap_end|gap_width|read_start|read_end|read_width|split_read_name|original_read_name"

    _name_re = re.compile(r"(.+)-([LR])-([\d]+)$")

    def __init__(self, split_read_name, chromosome, read_start, gap_start, 
            gap_end, read_end):
        self.chromosome = chromosome
        self.gap_start = gap_start
        self._gap_end = gap_end
        self._read_start = read_start
        self._read_end = read_end
        self._split_read_name = split_read_name
        self._hash_key = hash(tuple( \
            [chromosome, split_read_name, read_start, read_end]))

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False
        
    def __hash__(self):
        return self.hash_key
        
    def _read_width(self):
        return self._read_end - self._read_start
    
    def _gap_width(self):
        return self._gap_end - self.gap_start

    def _original_read_name(self):
        return self._name_re.match(self._split_read_name).group(1)

    def format(self, delimiter):
        return delimiter.join(
            [self.chromosome, str(self.gap_start), str(self._gap_end), 
                str(self._gap_width()), str(self._read_start), 
                str(self._read_end), str(self._read_width()), 
                self._split_read_name, self._original_read_name()])

class GapUtility():
    
    def __init__(self, original_read_len, delimiter):
        self._original_read_len = int(original_read_len)
        self._delimiter = delimiter

    def build_gap(self, sam_line):
        bits = sam_line.split(self._delimiter)[0:10]
        split_read_name = bits[0] 
        transcript_name = bits[2] 
        start_pos = bits[3]
        pnext = bits[7]
        seq = bits[9]
        l_gap_pos = int(start_pos) + len(seq)
        r_gap_pos = int(pnext)
        leftmost_start = int(start_pos)
        rightmost_end = int(pnext) + (self._original_read_len - len(seq))
        return Gap(split_read_name, transcript_name, leftmost_start, l_gap_pos,
            r_gap_pos, rightmost_end)

    def samfile_to_gaps(self, sam_file):
        gaps = []
        for line in sam_file:
            if line.startswith("@"):
                continue
            tlen = int(line.split(self._delimiter)[8])
            if tlen > 0: 
                gaps.append(self.build_gap(line))
        return gaps

    @staticmethod
    def sort_gaps(gaps):
        return sorted(gaps, key=lambda gap: (gap.chromosome, gap.gap_start))

    def write_gap_file(self, sorted_gaps, writer, additional_header_lines):
        for line in additional_header_lines:
            writer.write("#")
            writer.write(line)
            writer.write("\n")
        writer.write(Gap.header())
        writer.write("\n")
        for gap in sorted_gaps:
            writer.write(gap.format(self._delimiter))
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
    gaps = GapUtility.sort_gaps(gaps)

    logger.log("writing {0} gaps to file".format(len(gaps)))
    gap_file = open(gap_file_name, "w")
    gap_utility.write_gap_file(gaps, gap_file, header_lines)
    gap_file.close()
    
    
    logger.log("complete")
    
if __name__ == "__main__":
    BASENAME = os.path.basename(sys.argv[0])
    if (len(sys.argv) != 4):
        # pylint: disable=line-too-long
        print ("usage: {0} [sam_file] [original_read_len] [gap_file]".format(BASENAME))
        sys.exit()

    (SAM_FILE_NAME, ORIGINAL_READ_LEN, GAP_FILE_NAME) = sys.argv[1:]
    SAM_FILE_NAME = os.path.abspath(SAM_FILE_NAME)
    GAP_FILE_NAME = os.path.abspath(GAP_FILE_NAME)

    main(SAM_FILE_NAME, ORIGINAL_READ_LEN, GAP_FILE_NAME, "\t") 
    print ("{0} done.".format(BASENAME))
