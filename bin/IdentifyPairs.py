#! /usr/bin/env python

import sys, os
import datetime
import resource

class StdErrLogger():
	
	def __init__(self, verbose=0):
		self._verbose = verbose
	
	def log(self, message):
		print >> sys.stderr, "{0}|{1}".format(datetime.datetime.today(), message)

	
	def log_usage(self, message):
		self.log(message)
		if (self._verbose):
			usage = resource.getrusage(resource.RUSAGE_SELF)
			memory_used = (usage.ru_maxrss * resource.getpagesize())/(1024.0 * 1024)
			print >> sys.stderr, "{0}|usertime={1}|systime={2}|memory_used(mb)={3}".format(datetime.datetime.today(), usage.ru_utime, usage.ru_stime, memory_used)


class SplitRead():
	
	@classmethod
	def _swap_side(cls, side):
		if side == 'L':
			return 'R'
		else:
			return 'L'
	
	def __init__(self, read_len, name, side, split_len, strand, chr, position, matches):
		self._read_len = read_len
		self.name = name
		self.side = side
		self._split_len = split_len
		self._strand = strand
		self._chr = chr
		self._position = position
		self._matches = matches
	
	def _compliment_key(self):
		return (self.name, SplitRead._swap_side(self.side), (self._read_len - self._split_len), self._strand, self._chr)
	
	def key(self):
		if self.side == 'L':
			return (self.name, self.side, self._split_len, self._strand, self._chr)
		else:
			return self._compliment_key()
	
	def format(self, delimiter = "\t"):
		return delimiter.join([self.name, self.side, str(self._split_len), self._strand, self._chr, str(self._position), str(self._matches)])
	
	def distance(self, otherRead):
		return abs(self._position - otherRead._position)

class SplitReadBuilder():
	
	def __init__(self, read_len, delimiter = "\t"):
		self._read_len = read_len
		self._delimiter = delimiter
	
	def build(self, line):
		try:
			(name, side, split_len, strand, chr, position, seq, quality, matches) = line.split(self._delimiter)
			split_read = SplitRead(self._read_len, name, side, int(split_len), strand, chr, int(position), int(matches))
		except ValueError as e:
			raise SplitReadParseError(line, e)
		
		return split_read


class SplitReadParseError(Exception):
	def __init__(self, line, root_exception):
		self.line = line
		self.root_exception = root_exception
	
	def __str__(self):
		return repr("Could not parse line: '{0}' ; {1}".format(self.line, str(self.root_exception)))


def build_read_groups(split_read_builder, reader, logger):
	read_groups = {}
	count = 0
	for line in reader:
		count += 1
		if count % 100000 == 1:
			logger.log_usage("build_read_groups|processing line {0}".format(count))
		line = line.rstrip()
		split_read = split_read_builder.build(line)
		#if split_read._chr != "chr2":
		#	continue
		_add(read_groups, split_read)

	logger.log_usage("build_read_groups|done:processed {0} lines".format(count))
	return read_groups


def filter_orphans(read_groups, logger):
	count = 0
	for key in read_groups.keys():
		value = read_groups[key]
		if value[0] and value[1]:
			pass
		else:
			del read_groups[key]
		
		count += 1
		if count % 100000 == 1:
			logger.log_usage("prune_orphans|processing line {0}".format(count))
	
	logger.log("pruning orphan reads complete: {0} pairs remain".format(len(read_groups)))

def print_pairs(pairs, writer, delimiter, logger):
	count = 0
	for pair in pairs:
		count += 1
		left_read = pair[0].format(delimiter)
		right_read = pair[1].format(delimiter)
		distance = str(pair[2])
		writer.write(delimiter.join([left_read, right_read, distance]))
		writer.write("\n")
		if count % 100000 == 1:
			logger.log_usage("print_pairs|processing pair {0}".format(count))
	
	logger.log("print_pairs complete|processed {0} pairs".format(count))


def build_pairs_from_groups(read_groups, logger):
	count = 0
	pairs = []
	for read_group in read_groups.values():
		count += 1
		if count % 100000 == 1:
			logger.log_usage("build_pairs_from_groups|processing read_group {0}".format(count))
		left_reads = read_group[0]
		right_reads = read_group[1]
		for left_read in left_reads:
			for right_read in right_reads:
				distance = left_read.distance(right_read)
				pairs.append((left_read, right_read, distance))
	logger.log_usage("build_pairs_from_groups complete|processed {0} read_groups".format(count))
	return pairs

def filter_on_distance(pairs, min_distance, max_distance, logger):
	filtered_list = []
	for pair in pairs:
		distance = pair[2]
		if distance >= min_distance and distance <= max_distance:
		 filtered_list.append(pair)
	logger.log("filter_on_distance complete| {0} pairs processed, {1} pairs passed".format(len(pairs), len(filtered_list)))
	return filtered_list
	

def _add(d, split_read):
	key = split_read.key()
	if not key in d:
		d[key]=([],[])
	index = 0 if split_read.side == "L" else 1
	d[key][index].append(split_read)




def process_file(read_len, input_file_name, output_file_name, min_dist, max_dist):
	
	builder = SplitReadBuilder(read_len)
	reader = open(input_file_name, "r")
	writer = open(output_file_name, "w")
	logger = StdErrLogger()
	
	logger.log("process_file args| read_len:{0}, input_file_name:{1}, output_file_name:{2} minimum_distance:{3}, maximum_distance:{4}".format(read_len, \
				input_file_name, output_file_name, min_dist, max_dist))
	logger.log("processing file {0} begins".format(input_file_name))
	
	read_groups = build_read_groups(builder, reader, logger)
	filter_orphans(read_groups, logger)
	pairs = build_pairs_from_groups(read_groups, logger)
	pairs = filter_on_distance(pairs, min_dist, max_dist, logger)
	print_pairs(pairs, writer, "\t", logger)
	
	logger.log("processing file {0} complete".format(input_file_name))

if __name__ == "__main__":
	#TODO: adjust to use args/check_args

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

	process_file(read_len, infile, outfile, min_distance, max_distance) 
	print "done."
