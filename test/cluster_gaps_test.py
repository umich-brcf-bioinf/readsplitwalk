#! /usr/bin/env python2.6

import sys ; sys.path.insert(0, "../bin")
import unittest
from cluster_gaps import bed_line, samfile_to_bedfile


class ClusterGapsTest(unittest.TestCase):

	def test_bed_line(self):
		
		original_read_len = 50
		sam_line = "read_1-stuff-L-41|67|transcript42|100|score|cigar|=|150|50|ACGCT|CCCFF|XA:foo"
		field_delimiter = "|"
		name_delimiter = "#"
		actual_bed_line = bed_line(sam_line, original_read_len, field_delimiter, name_delimiter)
		self.assertEqual("transcript42|105|149|read_1-stuff-L-41#100#195\n", actual_bed_line)

	def test_samfile_to_bedfile(self):		
		original_read_len = 50
		sam_file = ["@header1", "@header2", 
				"read1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|stuff", 
				"read1|113|transcript42|200|score|cigar|=|150|-50|ACGCT|qual|stuff",
				"read2|113|transcript43|205|score|cigar|=|155|-50|ACGCT|qual|stuff",
				"read2|67|transcript43|155|score|cigar|=|205|50|ACGCT|qual|stuff"]

		writer = MockWriter()

		samfile_to_bedfile(sam_file, original_read_len, writer, "|")

		actual_lines = sorted(writer.lines())
		self.assertEqual(2, len(actual_lines))
		self.assertEqual(True, actual_lines[0].startswith("transcript42"))
		self.assertEqual(True, actual_lines[1].startswith("transcript43"))


class MockWriter():
	def __init__(self):
		self._content = []

	def write(self, content):
		self._content.append(content)
		
	def lines(self):
		return "".join(self._content).splitlines()

if __name__ == "__main__":
	unittest.main() 
