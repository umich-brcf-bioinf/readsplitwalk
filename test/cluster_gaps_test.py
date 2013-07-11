#! /usr/bin/env python2.6

import sys ; sys.path.insert(0, "../bin")
import unittest
from cluster_gaps import GapUtility


class ClusterGapsTest(unittest.TestCase):

	def test_gap_line(self):	
		original_read_len = 50
		delimiter = "|"
		sam_line = "read-name-L-42|67|transcript42|100|score|cigar|=|150|50|ACGCT|CCCFF|XA:foo"

		actual_bed_line = GapUtility(original_read_len, delimiter).gap_line(sam_line)
		
		actual_bits = actual_bed_line.split(delimiter)
		self.assertEqual(9, len(actual_bits), actual_bed_line)
		self.assertEqual("transcript42", actual_bits[0])
		self.assertEqual("105", actual_bits[1])
		self.assertEqual("150", actual_bits[2])
		self.assertEqual("45", actual_bits[3])
		self.assertEqual("100", actual_bits[4])
		self.assertEqual("195", actual_bits[5])
		self.assertEqual("95", actual_bits[6])
		self.assertEqual("read-name-L-42", actual_bits[7])
		self.assertEqual("read-name", actual_bits[8])

	def test_samfile_to_gap_lines_skipsNegativeStrand(self):		
		original_read_len = 50
		delimiter="|"
		sam_file = ["@header1", "@header2", 
				"read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|stuff", 
				"read1-L-1|113|transcript42|200|score|cigar|=|150|-50|ACGCT|qual|stuff",
				"read2-L-1|113|transcript43|205|score|cigar|=|155|-50|ACGCT|qual|stuff",
				"read2-L-1|67|transcript43|155|score|cigar|=|205|50|ACGCT|qual|stuff"]

		gap_lines = GapUtility(original_read_len, delimiter).samfile_to_gap_lines(sam_file)

		actual_lines = sorted(gap_lines)
		self.assertEqual(3, len(actual_lines))
		self.assertEqual("#chromosome|gapStart|gapEnd|gapWidth|readStart|readEnd|readWidth|splitReadName|originalReadName", actual_lines[0])
		self.assertEqual(True, actual_lines[1].startswith("transcript42"))
		self.assertEqual(True, actual_lines[2].startswith("transcript43"))

	def test_sort_gap_lines_sortedByChromByAlphaAndGapStartByNumeric(self):
		original_read_len = 50
		delimiter = "|"
		sorted_lines = ["chrom1|1|foo",
				"chrom10|1|bar",
				"chrom2|1|uno",
				"chrom2|2|dos",
				"chrom2|10|diez",
				"chrom2|11|once",
				"chrom20|1|froody"]

		input_lines = ["#header42"] 
		input_lines.extend(sorted_lines[::-1])
		input_lines.extend(["#header1"]) 

		actual_lines = GapUtility(original_read_len, delimiter).sort_gap_lines(input_lines)

		self.assertEqual(len(input_lines), len(actual_lines))
		self.assertEqual(["#header42","#header1"], actual_lines[:2])
		self.assertEqual(sorted_lines, actual_lines[2:])

	def test_write_gap_file(self):
		original_read_len = 50
		delimiter = "|"
		additional_header_lines=["hoopy", "frood"]
		input_lines = ["#foo", "bar","baz"]
		writer = MockWriter()

		GapUtility(original_read_len, delimiter).write_gap_file(input_lines, writer, additional_header_lines)

		actual_lines = writer.lines()
		expected_lines = ["#hoopy", "#frood"]
		expected_lines.extend(input_lines)
		self.assertEqual(expected_lines, actual_lines)


class MockWriter():
	def __init__(self):
		self._content = []

	def write(self, content):
		self._content.append(content)
		
	def lines(self):
		return "".join(self._content).splitlines()

if __name__ == "__main__":
	unittest.main() 
