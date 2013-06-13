#! /usr/bin/env python

#  This script will construct an executable shell script to process a set of files.

# ./IdentifyPairs.py <file> <outfile> <read_len> <min_dist> <max_dist> 2> <logfile> &

import sys, os, glob

if __name__ == "__main__":
	
	indir = sys.argv[1]
	scriptname = "./IdentifyPairs.py"
	filename_extension = sys.argv[2]
	read_len = 33
	min_dist = 2
	max_dist = 39999
	
	fullpath = os.path.abspath(indir)  
	if not os.path.isdir(indir):
		raise ValueError("indir [{0}] does not exist".format(fullpath))

	filespath = "{0}/*.{1}".format(fullpath, filename_extension)	
	filenames = glob.glob(filespath)
	filenames.sort()
	
	for fn in filenames:
		outfilename = fn+".pairs"
		logfilename = outfilename + ".log"
		cmdstring = "{0} {1} {2} {3} {4} {5} 2> {6} &".format(scriptname, fn, outfilename, read_len, min_dist, max_dist, logfilename)
		print cmdstring 