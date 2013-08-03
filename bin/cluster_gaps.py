#! /usr/bin/env python

"""
cluster_gaps.py
7/10/2013- cgates/pulintz
7/14/2013 - cgates: refactored to use a Gap object
7/28/2013 - cgates: overhauled to 
    a) use DBSCAN for clustering, 
    b) read and passthrough sample ids from sam 
    c) emit a sam file with cluster annotations
8/1/2013 - cgates: adjusted to emit original read as sam tag
"""
from contextlib import nested
import datetime
import os
import re
import resource
import sys
import traceback
from cluster_utility import DbscanClusterUtility

class ClusterGapsError(Exception):
    """Base class for exceptions in this module."""
    pass

class MissingReadGroupError(ClusterGapsError):
    def __init__(self, line):
        super(MissingReadGroupError, self).__init__()
        self.line = line

    def __str__(self):
        return repr("Alignment has no read group: '{0}'". \
            format(self.line))

class InvalidReadGroupError(ClusterGapsError):
    def __init__(self, line):
        super(InvalidReadGroupError, self).__init__()
        self.line = line

    def __str__(self):
        return repr("Specified read group not valid or not found in header: '{0}'". \
            format(self.line))

            
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
    def header(delimiter):
        return delimiter.join(["#chromosome", "cluster", "sample_name", 
            "gap_start", "gap_end","gap_width","read_start","read_end",
            "read_width","split_read_name","original_read_name"])

    _name_re = re.compile(r"(.+)-([LR])-([\d]+)$")

    def __init__(self, sample, split_read_name, chromosome, read_start, gap_start, 
            gap_end, read_end):
        self.sample = sample
        self.chromosome = chromosome
        self.gap_start = gap_start
        self._gap_end = gap_end
        self._read_start = read_start
        self._read_end = read_end
        self._split_read_name = split_read_name
        self.cluster = -1
        self._hash_key = hash(tuple( \
            [sample, chromosome, split_read_name, read_start, read_end]))

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False
        
    def __hash__(self):
        return self._hash_key
        
    def key(self):
        return self._hash_key
        
    def _read_width(self):
        return self._read_end - self._read_start
    
    def gap_width(self):
        return self._gap_end - self.gap_start

    def _original_read_name(self):
        return self._name_re.match(self._split_read_name).group(1)

    def format(self, delimiter):
        return delimiter.join(
            [self.chromosome, str(self.cluster), self.sample, str(self.gap_start), str(self._gap_end), 
                str(self.gap_width()), str(self._read_start), 
                str(self._read_end), str(self._read_width()), 
                self._split_read_name, self._original_read_name()])

    def additional_sam_tags(self, delimiter):
        return "XC:i:{0}{1}XR:Z:{2}".format(self.cluster, delimiter, self._original_read_name())        

class GapUtility():
    
    def __init__(self, original_read_len, delimiter, logger):
        self._original_read_len = int(original_read_len)
        self._delimiter = delimiter
        self._logger = logger
        self._read_group_sample_dict = {}

    def process_sam_header_line(self, sam_line):
        if sam_line.startswith("@RG"):
            key_values = sam_line.rstrip().split(self._delimiter)[1:]
            read_group = dict(u.split(":") for u in key_values)
            try:
                self._read_group_sample_dict[read_group['ID']] = read_group['SM']
            except:
                print (read_group)
                raise InvalidReadGroupError(sam_line)
                
    def sample_from_alignment(self, sam_line):
        bits = sam_line.rstrip().split(self._delimiter)
        if len(bits) < 12:
            raise MissingReadGroupError(sam_line)
        opt_dict = dict(re.split(':\w:', entries) for entries in bits[11:])
        if 'RG' not in opt_dict:
            raise MissingReadGroupError(sam_line)
        read_group = opt_dict['RG']
        if read_group not in self._read_group_sample_dict:
            raise InvalidReadGroupError(sam_line)
        return self._read_group_sample_dict[read_group]
        

    def build_gap(self, sam_line):
        bits = sam_line.split(self._delimiter)[0:10]
        split_read_name = bits[0] 
        transcript_name = bits[2] 
        start_pos = int(bits[3])
        pnext = int(bits[7])
        seq = bits[9]
        if start_pos < pnext :
            leftmost_start = start_pos
            gap_start = start_pos + len(seq)
            gap_end = pnext
            rightmost_end = pnext + (self._original_read_len - len(seq))
        else:
            leftmost_start = pnext
            gap_start = pnext + (self._original_read_len - len(seq))
            gap_end = start_pos
            rightmost_end = start_pos + len(seq)
        
        sample_name = self.sample_from_alignment(sam_line)
        
        return Gap(sample_name, split_read_name, transcript_name, 
            leftmost_start, gap_start, gap_end, rightmost_end)

    def samfile_to_gaps(self, sam_file):
        
        def is_leftmost_alignment_of_pair(line):
            tlen = int(line.split(self._delimiter)[8])
            return tlen > 0
        
        gaps = []
        for line in sam_file:
            if line.startswith("@"):
                self.process_sam_header_line(line)
            elif is_leftmost_alignment_of_pair(line):
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
        writer.write(Gap.header(self._delimiter))
        writer.write("\n")
        for gap in sorted_gaps:
            writer.write(gap.format(self._delimiter))
            writer.write("\n")

    def write_sam_file(
            self, input_sam_file, gaps, output_sam_file, 
            additional_header_lines):
            
        def _tagged_sam_line(line, gap_dict):
            sam_gap = self.build_gap(line)
            reference_gap = gap_dict[sam_gap.key()]
            additional_tags = reference_gap.additional_sam_tags(self._delimiter)
            return "{0}{1}{2}\n".format(
                    line.rstrip(), self._delimiter, additional_tags)
        
        self._logger.log("building gap dictionary")
        gap_dict = {}
        for gap in gaps:
            gap_dict[gap.key()] = gap
                
        for line in additional_header_lines:
                output_sam_file.write("@CO\t{0}\n".format(line))
        count = 0
        for line in input_sam_file:
            count += 1
            if count % 100000 == 1:
                self._logger.log("processing line {0}".format(count))
            if line.startswith("@"):
                output_sam_file.write(line)
            else:
                output_sam_file.write(_tagged_sam_line(line, gap_dict))
        self._logger.log("processed {0} lines".format(count))

def main(input_sam_file_name, original_read_len, gap_file_name, output_sam_file_name, delimiter):
    logger = StdErrLogger(verbose=True)
    logger.log(" ".join(sys.argv), verbose=False)
    header_lines = [str(datetime.datetime.today()), " ".join(sys.argv)] 
    gap_utility = GapUtility(original_read_len, delimiter, logger)
    
    logger.log("parsing sam file")
    with open(input_sam_file_name,"r") as sam_file:
        gaps = gap_utility.samfile_to_gaps(sam_file)

    logger.log("sorting {0} gaps".format(len(gaps)))    
    gaps = GapUtility.sort_gaps(gaps)

    logger.log("clustering gaps")
    cluster_utility = DbscanClusterUtility(logger=logger, deduplication_threshold=1000)
    cluster_utility.assign_clusters(gaps)

    logger.log("writing {0} gaps to file".format(len(gaps)))
    with open(gap_file_name, "w") as gap_file:
        gap_utility.write_gap_file(gaps, gap_file, header_lines)

    logger.log("writing sam file with clusters")
    with nested(open(input_sam_file_name,"r"), open(output_sam_file_name,"w")) \
            as (input_sam_file, output_sam_file):
        gap_utility.write_sam_file(
            input_sam_file, gaps, output_sam_file, header_lines)

    logger.log("{0} complete".format(input_sam_file_name))

if __name__ == "__main__":
    BASENAME = os.path.basename(sys.argv[0])
    if (len(sys.argv) != 5):
        # pylint: disable=line-too-long
        print ("usage: {0} [input_sam_file] [original_read_len] [gap_file] [output_sam_file]".format(BASENAME))
        sys.exit()

    (INPUT_SAM_FILE_NAME, ORIGINAL_READ_LEN, GAP_FILE_NAME, OUTPUT_SAM_FILE_NAME) = sys.argv[1:]
    INPUT_SAM_FILE_NAME = os.path.abspath(INPUT_SAM_FILE_NAME)
    GAP_FILE_NAME = os.path.abspath(GAP_FILE_NAME)
    OUTPUT_SAM_FILE_NAME = os.path.abspath(OUTPUT_SAM_FILE_NAME)

    main(INPUT_SAM_FILE_NAME, ORIGINAL_READ_LEN, GAP_FILE_NAME, OUTPUT_SAM_FILE_NAME, "\t") 
    print ("{0} done.".format(BASENAME))
