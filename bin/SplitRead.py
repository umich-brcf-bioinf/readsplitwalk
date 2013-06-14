#! /usr/bin/env python

import sys, os, re
import datetime
import resource


class FQStanza(object):
	""" Encapsulates the four line chunks of a fastq file. """
	
	@classmethod
	def parse(cls, line_str, delim="\n"):
		return FQStanza(*line_str.rstrip().split(delim))
		
	def __init__(self, main_header, seq, score_header, score):	
		(self.main_header, self.seq, self.score_header, self.score) = (main_header, seq, score_header, score)

	def split(self, split_position):
		right_size = len(self.seq)-split_position
		(left_seq, right_seq) = (self.seq[0:split_position], self.seq[split_position:])
		(left_score, right_score) = (self.score[0:split_position], self.score[split_position:])
		left_main_header = "{0}-L-{1}".format(self.main_header, split_position)
		left_score_header = "{0}-L-{1}".format(self.score_header, split_position)
		right_main_header = "{0}-R-{1}".format(self.main_header, right_size)
		right_score_header = "{0}-R-{1}".format(self.score_header, right_size)

		left_stanza = FQStanza(left_main_header, left_seq, left_score_header, left_score)
		right_stanza = FQStanza(right_main_header, right_seq, right_score_header, right_score)

		return (left_stanza, right_stanza)

	def __repr__(self):
		return "{0}\n{1}\n{2}\n{3}".format(self.main_header, self.seq, self.score_header, self.score)


def build_splits(in_stanza, split_margin):
	stanzas = []
	for split_position in range(split_margin, (len(in_stanza.seq)-split_margin)+1):
		stanzas.extend(in_stanza.split(split_position))
	return stanzas


def write_stanzas(in_stanzas, writer, split_margin):
	for in_stanza in in_stanzas:
		for out_stanza in build_splits(in_stanza, split_margin):
			writer.write(str(out_stanza))
			writer.write("\n")


def stanza_generator(reader, stanza_delimiter):
	#skip header
	for line in reader:
		stanza_str = line
		if stanza_str.startswith(stanza_delimiter):
			break
			
	for line in reader:
		if line.startswith(stanza_delimiter):
			yield FQStanza.parse(stanza_str)
			stanza_str = line
		else:
			stanza_str += line
	yield FQStanza.parse(stanza_str)

#def process_file():
#reader = open(file)
#generator = stanza_generator(reader)
#writer = open(file)
#write_stanzas(generator, kjdhsfk)

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