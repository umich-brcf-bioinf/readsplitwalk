#! /usr/bin/env python

import sys, os, re
import datetime
import resource


class FQStanza(object):
	""" Encapsulates the four line chunks of a fastq file. """
	
	def __init__(self, line_str, delim="\n"):
		(self.main_header, self.seq, self.score_header, self.score) = \
			line_str.rstrip().split(delim)
		
		



	
	
	
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