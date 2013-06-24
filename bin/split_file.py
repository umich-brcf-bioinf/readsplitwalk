#! /usr/bin/env python

import sys, os, re
import datetime
import resource

_SPLIT_FILE_EXTENSION = "prt"

def splitfile(filesystem, filepath, output_path, col_index, delim):
	count = 0
	# first pass to get values
	datafile = filesystem.open_file(filepath,"r")
	values = set()
	delim_re = re.compile(delim)
	for line in datafile:
		count += 1
		line = line.rstrip()
		value = delim_re.split(line)[col_index]
		values.add(value)

		if count % 10000 == 1:
			if len(values) > 100:
				raise ValueError("More than 100 distinct values for column index [{0}]".format(col_index))

	datafile.close()
	
	# construct subfiles
	subfiles = {}
	base_filename = os.path.basename(filepath)
	for value in values:
		subfilen = "{0}{1}.{2}.{3}".format(output_path, base_filename, value, _SPLIT_FILE_EXTENSION)
		subfile = filesystem.open_file(subfilen, "w")
		subfiles[value] = subfile
	
	# partition datafile
	datafile = filesystem.open_file(filepath, "r")
	for line in datafile:
		cvalue = delim_re.split(line.rstrip())[col_index]
		subfiles[cvalue].write(line)
	
	datafile.close()
	
	# cleanup
	for value in values:
		subfiles[value].close()
	
	
class FileSystem():
	
	def open_file(self, filename, mode):
		return open(filename, mode)
	
if __name__ == "__main__":

	if (len(sys.argv)  < 4 or len(sys.argv) > 5):
		print "usage: {0} [infile] [output path] [partition column zero-based index] [delimiter='\\t' (regex delimiter)]".format(os.path.basename(sys.argv[0]))
		sys.exit() 


	infile = sys.argv[1]
	if not os.path.isfile(infile):
		raise ValueError("infile [{0}] does not exist".format(infile))

	output_path = os.path.join(sys.argv[2],"") 
	if not os.path.isdir(output_path):
		raise ValueError("Output path [{0}] does not exist".format(output_path))

	col_index = int(sys.argv[3])

	if len(sys.argv) == 5:
		delim = sys.argv[4]
	else:
		 delim = "\t" 
	filesystem = FileSystem()
	splitfile(filesystem, infile, output_path, col_index, delim)
	print "done."