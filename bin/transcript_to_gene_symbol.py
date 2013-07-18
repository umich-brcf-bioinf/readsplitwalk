#! /usr/bin/env python
import sys, os

def parse_gene_map(infile_name):
	# transcipt -> gene_symbol
	infile = file(infile_name)
	gene_map = {}
	for line in infile.readlines():
		line = line.strip()
		bits = line.split("\t")
		if len(bits) < 3:
			print "No gene symbol for transcript: '{0}'".format(line)
			gene_map[bits[0]] = ""
			continue
		gene_map[bits[0]] = bits[2]
	infile.close()
	return gene_map
	

if __name__ == '__main__':
	
	if (len(sys.argv) != 5):
			print "usage: {0} [tsv transcript infile] [transcript_column_zero_based] [mapping datafile] [outfile]".format(os.path.basename(sys.argv[0]))
			sys.exit()
			
	transcript_filename = sys.argv[1]
	transcript_index=int(sys.argv[2])
	mapping_data = sys.argv[3]
	outfile_name = sys.argv[4]
	
	transcript_file = file(transcript_filename)
	outfile = file(outfile_name, "w")
	
	gene_map = parse_gene_map(mapping_data)
	
	missed_lines = 0
	total_transcript_lines = 0
	missed_keys = set()
	mapped_keys = set()
	for line in transcript_file.readlines():
		line = line.strip()
		total_transcript_lines += 1
		bits = line.split("\t")
		transcript_id = bits[transcript_index]
		gene_sym = "Not found."
		if transcript_id in gene_map:
			gene_sym = gene_map[transcript_id]
			mapped_keys.add(transcript_id)
		else:
			missed_lines += 1
			missed_keys.add(transcript_id)
			#print "No transcript mapping: {0}".format(transcript_id)
		outstr = "{0}\t{1}\n".format(line, gene_sym)
		outfile.write(outstr)
	outfile.close()
	print "{0} total transcript lines".format(total_transcript_lines)
	print "{0} mapped transcripts".format(len(mapped_keys))
	print "{0} missed transcript lines".format(missed_lines)
	print "{0} missed transcript IDs".format(len(missed_keys))
	print "done."
