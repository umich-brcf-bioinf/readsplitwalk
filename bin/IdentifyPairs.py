#! /usr/bin/env python

"""
IdentifyPairs.py
6/7/2013 - cgates/pulintz
Accepts a file of aligned split reads and creates a file of left/right read pairs based on a specified read length and min/max distance criteria. 
This ia a part of the read-split-walk process. 


Implementation:
The script parses the input file, assembling a set of "read groups keys". A "read group" is a set of matching left and right reads and the set of keys
define the universe of reads in the output (typically a small fraction of the input). 
Based on those keys, the program reads the file again, creating a hash of read groups which are then flattened to a collection of individual left-right pairs along with a distance (one line per pair). 
The collection of pairs is filtered on distance and written to the output file. 

This program is single threaded and will consume one processor for the duration of the run.
As of 6/16/2013, running a 16Gb input file (on ccmb-comp3 with no other load) took approximately 30 minutes and 30Gb of resident memory, producing a !Gb output file.  
Larger runs can be accomodated in less memory by partitioning the input files by chromosome. 


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
import sys


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
	def __init__(self, name, side, split_len, strand, chr, position, matches):
		self.name = name
		self.side = side
		self._split_len = split_len
		self._strand = strand
		self._chr = chr
		self._position = position
		self._matches = matches
	
	def format(self, delimiter = "\t"):
		return delimiter.join([self.name, self.side, str(self._split_len), self._strand, self._chr, str(self._position), str(self._matches)])
	
	def distance(self, otherRead):
		"""Absolute distance between this position and the otherRead"""
		return abs(self._position - otherRead._position)


class SplitReadBuilder():
	"""Interprets a key, side, or SplitRead from a line of text"""
	def __init__(self, read_len, delimiter = "\t"):
		self._read_len = read_len
		self._delimiter = delimiter
	
	def build(self, line):
		try:
			(name, side, split_len, strand, chr, position, seq, quality, matches) = line.rstrip().split(self._delimiter) 
			return SplitRead(name, side, int(split_len), strand, chr, int(position), int(matches))
		except ValueError as e:
			raise SplitReadParseError(line, e)
		
        def key_side(self, line, delimiter = "\t"):
		"""Parse a key and side form a line.  Note that "right side" reads return complementary "left side" keys."""
                (name, side, split_len, strand, chr) = line.split(self._delimiter)[:5]
                if side == "L":
                        return ("{0}|{1}|{2}|{3}|{4}".format(name, side, split_len, strand, chr), side)
                else:  
                        new_side = "L"
                        new_split_len = self._read_len - int(split_len)
                        return ("{0}|{1}|{2}|{3}|{4}".format(name, new_side, new_split_len, strand, chr), side)


class SplitReadParseError(Exception):
	def __init__(self, line, root_exception):
		self.line = line
		self.root_exception = root_exception
	
	def __str__(self):
		return repr("Could not parse line: '{0}' ; {1}".format(self.line, str(self.root_exception)))


def _identify_common_group_keys(split_read_builder, reader, logger):
	"""Reads every line, returning the set of all read keys that appeared on both the left and right sides.  
Each key in the result identifies a "read group"; all reads with these keys will appear in the output file"""

	#Circumvents a gc bug; see modifications.
	gc.disable()
	group_keys = { "L" : set(), "R" : set() }
	count = 0
	for line in reader:
		count +=1
		if count % 100000 == 1:
			logger.log("identify_group_keys|Processing line {0}".format(count))
		line = line.rstrip()
		(key, side) = split_read_builder.key_side(line)
		group_keys[side].add(key)
	logger.log("identify_group_keys|Processed {0} lines".format(count))
	
	logger.log("identify_group_keys|Intersecting {0} left keys with {1} right keys".format(len(group_keys["L"]), len(group_keys["R"]))) 
	common_keys = group_keys["L"].intersection(group_keys["R"])
	logger.log("identify_group_keys|Found {0} common keys".format(len(common_keys)))

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
		count += 1
		if count % 100000 == 1:
			logger.log("build_read_groups|processing line {0}".format(count))
		(key, side) = split_read_builder.key_side(line)
		if key in common_keys:
			_add(read_groups, key, split_read_builder.build(line))

	logger.log("build_read_groups|done:processed {0} lines".format(count))

	gc.enable()
	return read_groups


def _build_pairs_from_groups(read_groups, logger):
	"""For each group, generate the cartesian product of left and right pairs, returning a list of (left, right, distance) tuples."""
	gc.disable()
	count = 0
	pairs = []
	for read_group in read_groups.values():
		count += 1
		if count % 100000 == 1:
			logger.log("build_pairs_from_groups|processing read_group {0}".format(count))
		left_reads = read_group[0]
		right_reads = read_group[1]
		for left_read in left_reads:
			for right_read in right_reads:
				distance = left_read.distance(right_read)
				pairs.append((left_read, right_read, distance))
	logger.log("build_pairs_from_groups complete|processed {0} read_groups".format(count))
	gc.enable()
	return pairs


def _filter_on_distance(pairs, min_distance, max_distance, logger):
	filtered_list = []
	for pair in pairs:
		distance = pair[2]
		if distance >= min_distance and distance <= max_distance:
		 filtered_list.append(pair)
	logger.log("filter_on_distance complete| {0} pairs processed, {1} pairs passed".format(len(pairs), len(filtered_list)))
	return filtered_list


def _write_pairs(pairs, writer, logger, delimiter="\t"):
        count = 0
        for pair in pairs:
                count += 1
                left_read = pair[0].format(delimiter)
                right_read = pair[1].format(delimiter)
                distance = str(pair[2])
                writer.write(delimiter.join([left_read, right_read, distance]))
                writer.write("\n")
                if count % 100000 == 1:
                        logger.log("write_pairs|processing pair {0}".format(count))

        logger.log("write_pairs|processed {0} pairs".format(count))

	
def main(read_len, input_file_name, output_file_name, min_dist, max_dist):
	
	builder = SplitReadBuilder(read_len)
	logger = StdErrLogger(True)
	
	logger.log("process_file|read_len:{0}, input_file_name:{1}, output_file_name:{2} minimum_distance:{3}, maximum_distance:{4}".format(read_len, \
				input_file_name, output_file_name, min_dist, max_dist))
	logger.log("process_file|{0} begins".format(input_file_name))
	
	
	reader = open(input_file_name, "r")
	common_keys = _identify_common_group_keys(builder, reader, logger)
	reader.close()

	reader = open(input_file_name, "r")
	read_groups = _build_read_groups(common_keys, builder, reader, logger)
	reader.close()	

	pairs = _build_pairs_from_groups(read_groups, logger)
	pairs = _filter_on_distance(pairs, min_dist, max_dist, logger)

        writer = open(output_file_name, "w")	
	_write_pairs(pairs, writer, logger)
	writer.close()

	logger.log("process_file|output written to {0}".format(output_file_name))
	logger.log("process_file|{0} complete".format(input_file_name))

if __name__ == "__main__":

	if (len(sys.argv) != 6):
		print "usage: {0} [infile] [outfile] [read_len] [min_distance] [max_distance]".format(os.path.basename(sys.argv[0]))
		sys.exit() 

	(infile, outfile, read_len, min_distance, max_distance) = sys.argv[1:]

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

	main(read_len, infile, outfile, min_distance, max_distance) 
	print "done."