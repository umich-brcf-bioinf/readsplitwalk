#! /usr/bin/env python

"""
IdentifyPairs.py
6/7/2013 - cgates/pulintz
Accepts a file of aligned split reads and creates a file of left/right read pairs based on a specified read length and min/max distance criteria. 
This ia a part of the read-split-walk process. 

Example usage:
./IdentifyPairs.py alignedSplitReadsFromBowtie.txt pairedSplitReads.out 77 2 39999 2> pairedSplitReads.log
(See usage detals at bottom.)

Implementation:
The script parses the input file, assembling a set of "read groups keys". A "read group" is a set of matching left and right reads and the set of keys
define the universe of reads in the output (typically a small fraction of the input). 
Based on those keys, the program reads the file again, creating a hash of read groups which are then flattened to a collection of individual left-right pairs along with a distance (one line per pair). 
The collection of pairs is filtered on distance and written to the output file. 

This program is single threaded and will consume one processor for the duration of the run.
Running a 16Gb input file (on a large server with no other load) took approximately 30 minutes and 30Gb of resident memory, producing a 1Gb output file.  
Larger runs can be accomodated in less memory by partitioning the input files by chromosome.  See additional script SplitFile.py for safe ways to partition the input.


Modifications:
6/13/2013 - cgates
To reduce memory consumption and improve performance, added step to identify common read group keys prior to building split reads in memory. 

6/16/2013 - cgates
Disabled generational garbage collection within several methods to avoid a gc bug that causes appends to slow down geometrically in large lists.
Note this does not affect actual gc behavior as there are no cyclical data structures.
"""

import datetime
import gc
import os
import resource
import sys, re


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
			self._log("usertime(s)={0}|systime(s)={1}|peak_memory_used(mb)={2}|{3}".format(usage.ru_utime, usage.ru_stime, memory_used, message))
		else:
			self._log(message)


class SplitRead():
	"""Basic data structure for an individual read"""
	def __init__(self, name, side, split_len, strand, chr, position, matches, read_len, aligned=True):
		self.name = name
		self.side = side
		self.split_len = split_len
		self._strand = strand
		self._chr = chr
		self._position = position
		self._matches = matches
		self._read_len = read_len
		self.aligned = aligned
	
	def format(self, delimiter = "\t"):
		"""Returns a formatted string that accepts a field delimiter"""
		return delimiter.join([self.name, self.side, str(self.split_len), self._strand, self._chr, str(self._position), str(self._matches)])
		
	def distance(self, other):
		"""Absolute distance between this position and the otherRead"""
		return abs(self._position - other._position)

	def is_oriented(self, other):
		"""Returns true if other is on different side, same strand, and positioned correctly relative to self"""
		if self.side == other.side or self._strand != other._strand:
			return False
		(left_position, right_position)  = (self._position, other._position) if self.side == 'L' else (other._position, self._position)
		strand = 1 if self._strand == "+" else -1 
		return (right_position - left_position) * strand > 0

	def key(self):
		"""Returns a key for this read-group.  A left or right split read from the same initial read aligning to any position on the same chromosome and 
strand  would be part of the same read group.  For this reason, the key is always the 'left-side' key."""
		if self.side == "L":
			return "{0}|{1}|{2}|{3}|{4}".format(self.name, self.side, self.split_len, self._strand, self._chr)
		else:  
			new_side = "L"
			new_split_len = self._read_len - int(self.split_len)
			return "{0}|{1}|{2}|{3}|{4}".format(self.name, new_side, new_split_len, self._strand, self._chr)

	def __eq__(self, other):
		if type(other) is type(self):
			return self.__dict__ == other.__dict__
		return False


class LegacySplitReadBuilder():
	"""Interprets a SplitRead from a line of a non-standard text file."""
	def __init__(self, read_len, delimiter = "\t"):
		self._read_len = read_len
		self._delimiter = delimiter
	
	def build(self, line):
		try:
			(name, side, split_len, strand, chr, position, seq, quality, matches) = line.rstrip().split(self._delimiter)[:9]
			return SplitRead(name, side, int(split_len), strand, chr, int(position), int(matches), self._read_len)
		except ValueError as e:
			raise SplitReadParseError(line, e)

	def is_header(self, line):
		return False



class BowtieSplitReadBuilder():
	"""Interprets SplitRead from a line of a bowtie alignment file."""
	def __init__(self, read_len, delimiter = "\t"):
		self._read_len = read_len
		self._delimiter = delimiter
		self._name_re = re.compile(r"(.+)-([LR])-([\d]+)$")
		
	def build(self, line):
		try:
			(name, strand, chr, position, seq, quality, matches) = line.rstrip().split(self._delimiter)[:7]
			m = self._name_re.match(name)
			(subname, side, split_len) = (m.group(1), m.group(2), m.group(3))
			return SplitRead(subname, side, int(split_len), strand, chr, int(position), int(matches), self._read_len)
		except ValueError as e:
			raise SplitReadParseError(line, e)

	def is_header(self, line):
		return False


class SamSplitReadBuilder():
	"""Interprets SplitRead from a line of a SAM file."""
	def __init__(self, read_len, delimiter = "\t"):
		self._read_len = read_len
		self._delimiter = delimiter
		self._name_re = re.compile(r"(.+)-([LR])-([\d]+)$")
		
	def build(self, line):
		try:
			(name, flag, rname, position) = line.rstrip().split(self._delimiter)[:4]
			m = self._name_re.match(name)
			(subname, side, split_len) = (m.group(1), m.group(2), m.group(3))
			strand = "+" if int(flag) & 16 == 0 else "-"
			aligned = True if int(flag) & 4 == 0 else False
			return SplitRead(subname, side, int(split_len), strand, rname, int(position), None, self._read_len, aligned)
		except ValueError as e:
			raise SplitReadParseError(line, e)
			
	def is_header(self, line):
		return line.startswith("@")


def _identify_common_group_keys(split_read_builder, reader, logger, read_len):
	"""Reads every line, returning the set of all read keys that appeared on both the left and right sides.  
Each key in the result identifies a "read group"; the reads with these keys that pass other filtering criteria 
will appear in the output file"""

	#Circumvents a gc bug; see modifications.
	gc.disable()
	group_keys = { "L" : set(), "R" : set() }
	count = 0
	min_len = read_len
	max_len = 0
	for line in reader:
		#skip headers
		if split_read_builder.is_header(line): continue
		count +=1
		if count % 100000 == 1:
			logger.log("identify_group_keys|Processing line {0}".format(count))
		split_read = split_read_builder.build(line)
		if split_read.aligned:
			min_len = min(split_read.split_len, min_len)
			max_len = max(split_read.split_len, max_len)
			key = split_read.key()
			group_keys[split_read.side].add(key)
	logger.log("identify_group_keys|Processed {0} lines".format(count))
	
	logger.log("identify_group_keys|Intersecting {0} left keys with {1} right keys".format(len(group_keys["L"]), len(group_keys["R"]))) 
	common_keys = group_keys["L"].intersection(group_keys["R"])
	logger.log("identify_group_keys|Found {0} common keys".format(len(common_keys)))

	# test lengths
	computed_len = min_len + max_len 
	if (computed_len != read_len):
		raise ValueError("Read length validation failed.  Specified length {0} doesn't equal computed length {1}".format(read_len, computed_len))

	gc.enable()
	return common_keys


def _build_read_groups(common_keys, split_read_builder, reader, logger):
	"""Reads every line, returning a hash of read groups; each read group is a tuple of matching left and right reads. Only reads with key in common keys are included."""

	#Cirdumvents a gc bug; see modifications.
	gc.disable()

	def _add(d, key, split_read):
        	index = 0 if split_read.side == "L" else 1
        	group = d.setdefault(key, ([],[]))
        	group[index].append(split_read)

	read_groups = {}
	count = 0
	for line in reader:
		#skip headers
		if split_read_builder.is_header(line): continue
		count += 1
		if count % 100000 == 1:
			logger.log("build_read_groups|processing line {0}".format(count))
		split_read = split_read_builder.build(line)
		if split_read.aligned:
			key = split_read.key()
			if key in common_keys:
				_add(read_groups, key, split_read)
			
	logger.log("build_read_groups|done:processed {0} lines".format(count))

	gc.enable()
	return read_groups


def _build_pairs_from_groups(read_groups, logger):
	"""For each group, generate the cartesian product of left and right pairs, returning a list of (left, right, distance) tuples."""
	gc.disable()
	count = 0
	pairs = {}
	for key, read_group in read_groups.iteritems():
		pair = pairs.setdefault(key, [])
		count += 1
		if count % 100000 == 1:
			logger.log("build_pairs_from_groups|processing read_group {0}".format(count))
		left_reads = read_group[0]
		right_reads = read_group[1]
		for left_read in left_reads:
			for right_read in right_reads:
				distance = left_read.distance(right_read)
				pair.append((left_read, right_read, distance))
	logger.log("build_pairs_from_groups complete|processed {0} read_groups".format(count))
	gc.enable()
	return pairs


def _filter_on_distance(all_read_group_pairs, min_distance, max_distance, logger):
	filtered_pairs = {}
	count_total = 0
	count_included = 0
	for key, pairs in all_read_group_pairs.iteritems():
		pair_list = []
		for pair in pairs:
			count_total += 1
			distance = pair[2]
			if distance >= min_distance and distance <= max_distance:
		 		pair_list.append(pair)
				count_included += 1
		if pair_list:
			filtered_pairs[key]=pair_list

	logger.log("filter_on_distance complete| {0} pairs processed, {1} pairs passed".format(count_total, count_included))
	return filtered_pairs


def _filter_on_orientation(all_read_group_pairs, logger):
        filtered_pairs = {}
        count_total = 0
        count_included = 0
        for key, pairs in all_read_group_pairs.iteritems():
                pair_list = []
                for pair in pairs:
                        count_total += 1
                        if pair[0].is_oriented(pair[1]):
                                pair_list.append(pair)
                                count_included += 1
                if pair_list:
                        filtered_pairs[key]=pair_list

        logger.log("filter_on_orientation complete| {0} pairs processed, {1} pairs passed".format(count_total, count_included))
        return filtered_pairs


def _write_rsw_pairs(all_read_group_pairs, writer, logger, delimiter="\t"):
        count = 0
        for read_group_pairs in all_read_group_pairs.itervalues():
		for pair in read_group_pairs:
                	count += 1
                	left_read = pair[0].format(delimiter)
                	right_read = pair[1].format(delimiter)
                	distance = str(pair[2])
                	writer.write(delimiter.join([left_read, right_read, distance]))
                	writer.write("\n")
                	if count % 100000 == 1:
                        	logger.log("write_rsw_pairs|processing pair {0}".format(count))

        logger.log("write_rsw_pairs|processed {0} pairs".format(count))


def _write_sam_pairs(read_group_pairs, reader, builder, writer, logger):
	count = 0
	for line in reader:
		count += 1
		if builder.is_header(line):
			writer.write(line)
		else:
			split_read = builder.build(line)
			if split_read.aligned and split_read.key() in read_group_pairs:
				for pair in read_group_pairs[split_read.key()]:
					if split_read == pair[0] or split_read == pair[1]:
						writer.write(line)

		if count % 100000 == 1:
				logger.log("write_sam_pairs|processing line {0}".format(count))

	logger.log("write_sam_pairs|processed {0} lines".format(count))
	
	
def main(read_len, input_file_name, output_file_name, sam_output_file_name, min_dist, max_dist):
	
	builder = SamSplitReadBuilder(read_len)
	logger = StdErrLogger(True)
	
	logger.log("process_file|read_len:{0}, input_file_name:{1}, output_file_name:{2}, sam_output_file_name:{3}, minimum_distance:{4}, maximum_distance:{5}".format(read_len, \
				input_file_name, output_file_name, sam_output_file_name, min_dist, max_dist))
	logger.log("process_file|{0} begins".format(input_file_name))
	
	
	reader = open(input_file_name, "r")
	common_keys = _identify_common_group_keys(builder, reader, logger, read_len)
	reader.close()

	reader = open(input_file_name, "r")
	read_groups = _build_read_groups(common_keys, builder, reader, logger)
	reader.close()	

	read_group_pairs = _build_pairs_from_groups(read_groups, logger)
	read_group_pairs = _filter_on_distance(read_group_pairs, min_dist, max_dist, logger)
	read_group_pairs = _filter_on_orientation(read_group_pairs, logger)

	writer = open(output_file_name, "w")	
	_write_rsw_pairs(read_group_pairs, writer, logger)
	writer.close()

	reader = open(input_file_name, "r")	
	writer = open(sam_output_file_name, "w")	
	_write_sam_pairs(read_group_pairs, reader, builder, writer, logger)
	writer.close()
	reader.close()

	logger.log("process_file|output written to {0}".format(output_file_name))
	logger.log("process_file|{0} complete".format(input_file_name))
	

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
