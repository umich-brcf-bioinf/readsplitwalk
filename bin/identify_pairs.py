#! /usr/bin/env python

"""
identify_pairs.py
6/7/2013 - cgates/pulintz
Accepts a file of aligned split reads and creates two output files with left/right read pairs based on a specified read length and min/max distance criteria. 
The first output file is a custom format which shows each pair on a single line including a distance between the two alignments.
The second output file is a subset of the input file (each alignment on a single line) filtered to contain only paired alignments. 


Example usage:
./identify_pairs.py alignedSplitReadsFromBowtie.sam pairedSplitReads.out 77 2 39999 2> pairedSplitReads.log
(See usage detals at bottom.)


Implementation:
The script parses the input file, assembling a set of "read groups keys". A "read group" is a set of matching left and right reads and the set of keys
define the universe of reads in the output (typically a small fraction of the input). 
Based on those keys, the program reads the file again, creating a hash of read groups which are then flattened to a collection of individual left-right pairs along with a distance (one line per pair). 
The collection of pairs is filtered on distance and written to the output file. 

This program is single threaded and will consume one processor for the duration of the run.
Running a 50Gb SAM input file (on a large server with no other load) took approximately 3 hours and 15Gb of resident memory, producing a 200Mb output file.  

The input file can be all the alignments for a sample. Processing can be paralellized (and memory footprint reduced) by splitting input into multiple files,
partitioning by chromosome, strand, or both. Importantly, the program assumes that all possible pairs exist within a single file, so each file should contain 
all alignments for a chromosome/strand. Thus, for example, you cannot partition the input files by chromosome region.
See partition_file.py for details on partitioning a file.

Modifications:
6/13/2013 - cgates
To reduce memory consumption and improve performance, added step to identify common read group keys prior to building split reads in memory. 

6/16/2013 - cgates
Disabled generational garbage collection within several methods to avoid a gc bug that causes appends to slow down geometrically in large lists.
Note this does not affect actual gc behavior as there are no cyclical data structures.

6/21/2013 - cgates/pulintz
Added ability to accept SAM file as input

6/24/2013 - cgates/pulintz
Added generation of SAM file alongside original output file. 

6/26/2013 - cgates
Adjusted SplitRead to not store sequence/quality (to reduce memory usage).
Adjusted to skip non-aligned reads in a SAM input file.
Adjusted generation of filtered SAM file to pass through subset of input lines instead of creating SAM data from SplitRead objects.
Added filter on pair orientation to exclude pairs whose left/right-strand orientation was incorrect.


6/27/2013 - cgates/pulintz
Switched from if-blocks to polymorphic handling of unaligned reads by introducing two new classes UnalignedSplitRead, SplitReadValidator, and also 
refactoring methods into SplitRead (write_sam_pairs, add_to_read_groups, add_to_group_keys, check_split_length). This reduced method sizes, improved SplitRead encapsulation,
simplified testing, and (because we now avoid creating a SplitRead for unaligned reads), improved performance.
Refatored how filtering works to use a filter chain/composite filter pattern.
Minor cosmetic adjustments to Logger.

7/01/2013 - pulintz
Added logic to SplitRead.write_sam_pairs to export split reads as paired reads, which permits IGV to use its paired-read display functionality.  Added an attribute
to SplitRead called _original_name to facilitate this.

7/8/2013 - cgates/pulintz
Renamed SplitRead.distance() to gap_distance and revised to calculate gap distance correctly.
"""

import datetime
import gc
import os
import re
import resource	
import sys
import traceback


class IdentifyPairsException(Exception):
	"""Base class for exceptions in this module."""
	pass

class SplitReadParseError(IdentifyPairsException):
	def __init__(self, line, root_exception):
		self.line = line
		self.root_exception = root_exception
	
	def __str__(self):
		return repr("Could not parse line: '{0}' ; {1}".format(self.line, str(self.root_exception)))


class StdErrLogger():
	"""Writes basic utilization data to stderr"""	
	def __init__(self, verbose = False):
		self._verbose = verbose
	
	def _log(self, message):
		print >> sys.stderr, "{0}|{1}".format(datetime.datetime.today(), message)

	
	def log(self, message):
		if (self._verbose):
			usage = resource.getrusage(resource.RUSAGE_SELF)
			memory_used = usage.ru_maxrss/1024
			function_name = traceback.extract_stack()[-2:-1][0][2]	
			self._log("usertime(s)={0:.0f}|systime(s)={1:.0f}|peak_memory_used(mb)={2}|{3}|{4}".format(usage.ru_utime, usage.ru_stime, memory_used, function_name, message))
		else:
			self._log(message)


class SplitRead():
	"""Basic data structure for an individual read"""
	def __init__(self, name, side, split_len, strand, chr, position, matches, original_read_len):
		self._name = name
		self._side = side
		self._split_len = split_len
		self._strand = strand
		self._chr = chr
		self._position = position
		self._matches = matches
		self._original_read_len = original_read_len

	def format(self, delimiter = "\t"):
		"""Returns a formatted string that accepts a field delimiter"""
		return delimiter.join([self._name, self._side, str(self._split_len), self._strand, self._chr, str(self._position), str(self._matches)])

	def gap_distance(self, other):
		"""Returns distance between the end of one read and the beginning of the other.  (Assumes all alignment positions are leftmost coordinate of positive strand.)"""
		(leftmost_read, rightmost_read) = (self, other) if self._position < other._position else (other, self)
		return rightmost_read._position - (leftmost_read._position + leftmost_read._split_len)

	def is_oriented(self, other):
		"""Returns true if other is on different side, same strand, and positioned correctly relative to self"""
		if self._side == other._side or self._strand != other._strand:
			return False
		(left_position, right_position)  = (self._position, other._position) if self._side == 'L' else (other._position, self._position)
		strand = 1 if self._strand == "+" else -1 
		return (right_position - left_position) * strand > 0

	def key(self):
		"""Returns a key for this read-group.  A left or right split read from the same initial read aligning to any position on the same chromosome and 
strand  would be part of the same read group.  For this reason, the key is always the 'left-side' key."""
		if self._side == "L":
			return "{0}|{1}|{2}|{3}|{4}".format(self._name, self._side, self._split_len, self._strand, self._chr)
		else:  
			new_side = "L"
			new_split_len = self._original_read_len - int(self._split_len)
			return "{0}|{1}|{2}|{3}|{4}".format(self._name, new_side, new_split_len, self._strand, self._chr)

	def left_name(self):
		"""Returns the left handed name."""
		if self._side == "L":
			new_split_len = self._split_len
		else:
			new_split_len = self._original_read_len - int(self._split_len)
		
		return "{0}-{1}-{2}".format(self._name, "L", new_split_len)

	def __eq__(self, other):
		if type(other) is type(self):
			return self.__dict__ == other.__dict__
		return False

	def write_sam_pairs(self, read_group_pairs, line, writer, delimiter = "\t"):
		key = self.key()
		if not key in read_group_pairs:
			return
		bits = line.split(delimiter)
		(chrom, mapq, cigar, pos, seq, qual, extra) = (bits[2], bits[4], bits[5], bits[7], bits[9], bits[10], bits[11:])
		extras = " ".join(extra)
		flag = pnext = tlen = 0
		for pair in read_group_pairs[key]:
			if self == pair[0]:
				pnext = pair[1]._position
				tlen = pnext-self._position
				if self._position < pnext:
					flag = 67
				else:
					flag = 115
			elif self == pair[1]:
				pnext = pair[0]._position
				tlen = pnext-self._position
				if self._position > pnext:
					flag = 131
				else:
					flag = 147
			else:
				continue	
			outstring = delimiter.join([self.left_name(), str(flag), chrom, str(self._position), mapq, cigar, "=", str(pnext), str(tlen), seq, qual, extras])
			writer.write(outstring)
				

	def add_to_read_groups(self, common_keys, read_groups):	
		key = self.key()
		if key in common_keys:
			index = 0 if self._side == "L" else 1
			group = read_groups.setdefault(key, ([],[]))
			group[index].append(self)

	def add_to_group_keys(self, group_keys):
		group_keys[self._side].add(self.key())

	def check_split_length(self, validator):
		validator.check_split_length(self._split_len)


class UnalignedSplitRead():
	"""No-op split read. Builders return a singleton instance of this class when buidling a read which has not aligned."""

	def format(self, delimiter=""): return None
	
	def distance(self, other): return None

	def is_oriented(self, other): return None

	def key(self): return None
		
	def write_sam_pairs(self, read_group_pairs, line, writer): pass

	def add_to_read_groups(self, common_keys, read_groups): pass
		
	def add_to_group_keys(self, group_keys): pass

	def check_split_length(self, validator): pass


class ReadLengthValidationError(IdentifyPairsException):
	pass

class ReadLengthValidator():
	"""Validates that no split_length exceeds the read_length and that the min + max split_lengths match the read_length.  Raises error if anything is amiss."""
        def __init__(self, original_read_len):
		self._original_read_len = original_read_len
                self._min_len = 1000000
                self._max_len = 0

        def check_split_length(self, split_len):
		if split_len > self._original_read_len:
			raise ReadLengthValidationError("Length of a split read ({0}) exceeds specified overall read length ({1})".format(split_len, self._original_read_len))
                self._min_len = min(split_len, self._min_len)
                self._max_len = max(split_len, self._max_len)

        def check_read_length(self):
                computed_len = self._min_len + self._max_len
                if (computed_len != self._original_read_len):
                        raise ReadLengthValidationError("Specified read length ({0}) doesn't equal computed read length ({1}) ({2}-{3})".format(self._original_read_len, computed_len, self._min_len, self._max_len))

	
class LegacySplitReadBuilder():
	"""Interprets a SplitRead from a line of a non-standard text file."""
	def __init__(self, original_read_len, delimiter = "\t"):
		self._original_read_len = original_read_len
		self._delimiter = delimiter
	
	def build(self, line):
		try:
			(name, side, split_len, strand, chr, position, seq, quality, matches) = line.rstrip().split(self._delimiter)[:9]
			return SplitRead(name, side, int(split_len), strand, chr, int(position), int(matches), self._original_read_len)
		except ValueError as e:
			raise SplitReadParseError(line, e)

	def is_header(self, line):
		return False


class BowtieSplitReadBuilder():
	"""Interprets SplitRead from a line of a bowtie alignment file."""
	def __init__(self, original_read_len, delimiter = "\t"):
		self._original_read_len = original_read_len
		self._delimiter = delimiter
		self._name_re = re.compile(r"(.+)-([LR])-([\d]+)$")
		
	def build(self, line):
		try:
			(name, strand, chr, position, seq, quality, matches) = line.rstrip().split(self._delimiter)[:7]
			m = self._name_re.match(name)
			(subname, side, split_len) = (m.group(1), m.group(2), m.group(3))
			return SplitRead(subname, side, int(split_len), strand, chr, int(position), int(matches), self._original_read_len)
		except ValueError as e:
			raise SplitReadParseError(line, e)

	def is_header(self, line):
		return False


class SamSplitReadBuilder():
	"""Interprets SplitRead from a line of a SAM file."""
	def __init__(self, original_read_len, delimiter = "\t"):
		self._original_read_len = original_read_len
		self._delimiter = delimiter
		self._name_re = re.compile(r"(.+)-([LR])-([\d]+)$")
		self.UNALIGNED_READ = UnalignedSplitRead()
		
	def build(self, line):
		try:
			(name, flag, rname, position) = line.rstrip().split(self._delimiter)[:4]
			aligned = int(flag) & 4 == 0
			
			if not aligned:
				 return self.UNALIGNED_READ

			m = self._name_re.match(name)
			(subname, side, split_len) = (m.group(1), m.group(2), m.group(3))
			strand = "+" if int(flag) & 16 == 0 else "-"
			return SplitRead(subname, side, int(split_len), strand, rname, int(position), None, self._original_read_len)

		except ValueError as e:
			raise SplitReadParseError(line, e)
			
	def is_header(self, line):
		return line.startswith("@")


def _identify_common_group_keys(split_read_builder, validator, reader, logger):
	"""Reads every line, returning the set of all read keys that appeared on both the left and right sides.  
Each key in the result identifies a "read group"; the reads with these keys that pass other filtering criteria 
will appear in the output file"""

	#Circumvents a gc bug; see modifications.
	gc.disable()
	group_keys = { "L" : set(), "R" : set() }
	count = 0

	for line in reader:
		if split_read_builder.is_header(line): continue
		count +=1
		if count % 100000 == 1: logger.log("processing line {0}".format(count))
		split_read = split_read_builder.build(line)
		split_read.add_to_group_keys(group_keys)
		split_read.check_split_length(validator)

	logger.log("processed {0} lines".format(count))
	
	logger.log("intersecting {0} left keys with {1} right keys".format(len(group_keys["L"]), len(group_keys["R"]))) 
	common_keys = group_keys["L"].intersection(group_keys["R"])
	logger.log("found {0} common keys".format(len(common_keys)))

	validator.check_read_length()

	gc.enable()
	return common_keys


def _build_read_groups(common_keys, split_read_builder, reader, logger):
	"""Reads every line, returning a hash of read groups; each read group is a tuple of matching left and right reads. Only reads with key in common keys are included."""

	#Circumvents a gc bug; see modifications.
	gc.disable()

	read_groups = {}
	count = 0
	for line in reader:
		if split_read_builder.is_header(line): continue
		count += 1
		if count % 100000 == 1: logger.log("processing line {0}".format(count))
		split_read = split_read_builder.build(line)
		split_read.add_to_read_groups(common_keys, read_groups)

	logger.log("processed {0} lines".format(count))

	gc.enable()
	return read_groups


def _build_pairs_from_groups(read_groups, logger):
	"""For each group, generate the cartesian product of left and right pairs, returning a dict of split_read_key, list of (left, right, distance) tuples."""
	gc.disable()
	count = 0
	pairs = {}
	for key, read_group in read_groups.iteritems():
		pair = pairs.setdefault(key, [])
		count += 1
		if count % 100000 == 1: logger.log("processing read_group {0}".format(count))
		left_reads = read_group[0]
		right_reads = read_group[1]
		for left_read in left_reads:
			for right_read in right_reads:
				pair.append((left_read, right_read))
	logger.log("processed {0} read_groups".format(count))
	gc.enable()
	return pairs


def _filter_pairs(all_read_group_pairs, filter, logger):
	"""Iterates over all pairs in the read group dict applying the specified filter and returning a new dict of the filtered results."""
	filtered_pairs = {}
	count_total = 0
	count_included = 0
	for key, pairs in all_read_group_pairs.iteritems():
		pair_list = []
		for pair in pairs:
			count_total += 1
			if filter(pair[0], pair[1]):				
		 		pair_list.append(pair)
				count_included += 1
		if pair_list:
			filtered_pairs[key]=pair_list

	logger.log("{0} pairs processed, {1} pairs passed".format(count_total, count_included))
	return filtered_pairs


def _distance_filter(min_distance, max_distance):
	def filter_pair(read1, read2):
		distance = read1.gap_distance(read2)
                return distance >= min_distance and distance <= max_distance;
	return filter_pair

def _orientation_filter(read1, read2):
	return read1.is_oriented(read2)

def _composite_filter(filter_list):
	def filter_pair(read1, read2):
		include = True
		for filter in filter_list:
			include = include and filter(read1, read2)
		return include
	return filter_pair

def _write_rsw_pairs(all_read_group_pairs, writer, logger, delimiter="\t"):
        count = 0
        for read_group_pairs in all_read_group_pairs.itervalues():
		for pair in read_group_pairs:
                	count += 1
                	left_read = pair[0].format(delimiter)
                	right_read = pair[1].format(delimiter)
                	distance = str(pair[0].gap_distance(pair[1]))
                	writer.write(delimiter.join([left_read, right_read, distance]))
                	writer.write("\n")
                	if count % 100000 == 1:
                        	logger.log("processing pair {0}".format(count))

        logger.log("processed {0} pairs".format(count))


def _write_sam_pairs(read_group_pairs, reader, builder, writer, logger, delimiter="\t"):
	count = 0
	for line in reader:
		count += 1
		if builder.is_header(line):
			writer.write(line)
		else:
			split_read = builder.build(line)
			split_read.write_sam_pairs(read_group_pairs, line, writer)

		if count % 100000 == 1:
				logger.log("processing line {0}".format(count))

	logger.log("processed {0} lines".format(count))
	
	
def main(original_read_len, input_file_name, output_file_name, sam_output_file_name, min_dist, max_dist):
	
	logger = StdErrLogger(True)
	logger.log("read_len:{0}, input_file_name:{1}, output_file_name:{2}, sam_output_file_name:{3}, minimum_distance:{4}, maximum_distance:{5}".format(read_len, \
				input_file_name, output_file_name, sam_output_file_name, min_dist, max_dist))
	logger.log("{0} begins".format(input_file_name))
	
	builder = SamSplitReadBuilder(original_read_len)
	validator = ReadLengthValidator(original_read_len)
	
	reader = open(input_file_name, "r")
	common_keys = _identify_common_group_keys(builder, validator, reader, logger)
	reader.close()

	reader = open(input_file_name, "r")
	read_groups = _build_read_groups(common_keys, builder, reader, logger)
	reader.close()	

	read_group_pairs = _build_pairs_from_groups(read_groups, logger)

	filter = _composite_filter([_distance_filter(min_dist, max_dist), _orientation_filter])
	read_group_pairs = _filter_pairs(read_group_pairs, filter, logger)
	 
	writer = open(output_file_name, "w")	
	_write_rsw_pairs(read_group_pairs, writer, logger)
	writer.close()

	reader = open(input_file_name, "r")	
	writer = open(sam_output_file_name, "w")	
	_write_sam_pairs(read_group_pairs, reader, builder, writer, logger)
	writer.close()
	reader.close()

	logger.log("output written to {0}".format(output_file_name))
	logger.log("{0} complete".format(input_file_name))
	

if __name__ == "__main__":

	if (len(sys.argv) != 6):
		print "usage: {0} [infile] [outfile] [read_len] [min_distance] [max_distance]".format(os.path.basename(sys.argv[0]))
		sys.exit() 

	(infile, outfile, read_len, min_distance, max_distance) = sys.argv[1:]
	infile = os.path.abspath(infile)
	outfile = os.path.abspath(outfile)
	sam_outfile = "{0}.sam".format(os.path.splitext(outfile)[0])

	# check params
	try:
		min_distance = int(min_distance)
		max_distance = int(max_distance)
		read_len = int(read_len)
		if (max_distance - min_distance < 1):
			raise ValueError("max distance must be greater than min distance")
	except ValueError as e:
		print str(e)
		print "usage: {0} [infile] [outfile] [read_len] [min_distance] [max_distance]".format(os.path.basename(sys.argv[0]))
		sys.exit() 	

	main(read_len, infile, outfile, sam_outfile, min_distance, max_distance) 
	print "done."
