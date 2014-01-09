#! /usr/bin/env python
""" Quick script to add readgroup lines to sam files. Script also appends the optional 'RG' tag data to the end of each data line. """

import sys, os
    
if __name__ == '__main__':
    
	if (len(sys.argv) != 5):
		print "usage: {0} [sam infile] [ReadGroup ID] [ReadGroup sample] [outfile]".format(os.path.basename(sys.argv[0]))
		sys.exit()
            
	sam_filename = sys.argv[1]
	rg_id = sys.argv[2]
	rg_sample = sys.argv[3]
	outfile_name = sys.argv[4]

	newline = "@RG\tID:{0}\tSM:{1}\n".format(rg_id, rg_sample)
	
	sam_file = file(sam_filename)
	outfile = file(outfile_name, "w")

	for line in sam_file.readlines():
		if line.startswith("@PG"):
			outfile.write(newline)
			outfile.write(line)
		elif line.startswith("@"):
			outfile.write(line)
		else:
			line = line.strip()
			newline = "{0}\tRG:Z:{1}\n".format(line, rg_id)
			outfile.write(newline)
    
	print "done."
