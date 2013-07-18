import unittest
from bin.cluster_gaps import GapUtility, Gap


class GapTestCase(unittest.TestCase):

	def test_gap_width(self):
		self.assertEqual(12, Gap("split-read-name-L-13", "chromosome", 0, 4, 16, 64)._gap_width())

	def test_read_width(self):
		self.assertEqual(64, Gap("split-read-name-L-13", "chromosome", 0, 4, 16, 64)._read_width())

	def test_original_read_name(self):
		self.assertEqual("split-read-name", Gap("split-read-name-L-13", "chromosome", 0, 4, 16, 64)._original_read_name())

	def test_header(self):
		self.assertEqual("#chromosome|gap_start|gap_end|gap_width|read_start|read_end|read_width|split_read_name|original_read_name", Gap.header())
		
	def test_format(self):
		gap = Gap("split-read-name-L-13", "chromosome", 0, 4, 16, 64)
		self.assertEqual("chromosome|4|16|12|0|64|64|split-read-name-L-13|split-read-name", gap._format("|"))

	def test_eq(self):
		base = Gap("split-read-name-L-13", "chromosome", 0, 4, 16, 64)
		self.assertEqual(True, base == Gap("split-read-name-L-13", "chromosome", 0, 4, 16, 64))
		self.assertEqual(False, base == "foo")
		self.assertEqual(False, base == Gap("X", "chromosome", 0, 4, 16, 64))
		self.assertEqual(False, base == Gap("split-read-name-L-13", "chromosome", -42, 4, 16, 64))
		self.assertEqual(False, base == Gap("split-read-name-L-13", "chromosome", 0, -42, 16, 64))
		self.assertEqual(False, base == Gap("split-read-name-L-13", "chromosome", 0, 4, -42, 64))
		self.assertEqual(False, base == Gap("split-read-name-L-13", "chromosome", 0, 4, 16, -42))


class GapUtilityTestCase(unittest.TestCase):

	def test_build_gap(self):	
		original_read_len = 50
		delimiter = "|"
		sam_line = "read-name-L-42|67|transcript42|100|score|cigar|=|150|50|ACGCT|CCCFF|XA:foo"

		gap = GapUtility(original_read_len, delimiter).build_gap(sam_line)
		
		self.assertEqual("transcript42", gap._chromosome)
		self.assertEqual(105, gap._gap_start)
		self.assertEqual(150, gap._gap_end)
		self.assertEqual(100, gap._read_start)
		self.assertEqual(195, gap._read_end)
		self.assertEqual("read-name-L-42", gap._split_read_name)

	def test_samfile_to_gaps_skipsNegativeStrand(self):		
		original_read_len = 50
		delimiter="|"
		sam_file = ["@header1", "@header2", 
				"read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|stuff", 
				"read1-L-1|113|transcript42|200|score|cigar|=|150|-50|ACGCT|qual|stuff",
				"read2-L-1|113|transcript43|205|score|cigar|=|155|-50|ACGCT|qual|stuff",
				"read2-L-1|67|transcript43|155|score|cigar|=|205|50|ACGCT|qual|stuff"]

		gaps = GapUtility(original_read_len, delimiter).samfile_to_gaps(sam_file)

		self.assertEqual(2, len(gaps))
		self.assertEqual("transcript42", gaps[0]._chromosome)
		self.assertEqual("transcript43", gaps[1]._chromosome)

	def test_sort_gap_lines_sortedByChromByAlphaAndGapStartByNumeric(self):
		
		original_read_len = 50
		delimiter = "|"
		sorted_gaps = [
			init_gap("chrom1", 1, "foo"),
			init_gap("chrom10", 1, "bar"),
			init_gap("chrom2", 1, "uno"),
			init_gap("chrom2", 2, "dos"),
			init_gap("chrom2", 10, "diez"),
			init_gap("chrom2", 11, "once"),
			init_gap("chrom20", 1, "froody")]

		input_gaps = sorted_gaps[::-1]

		actual_gaps = GapUtility(original_read_len, delimiter).sort_gaps(input_gaps)

		self.assertEqual(len(input_gaps), len(actual_gaps))
		self.assertEqual(sorted_gaps, actual_gaps)

	def test_write_gap_file(self):
		original_read_len = 50
		delimiter = "|"
		additional_header_lines=["hoopy", "frood"]
		gaps = [MockGap("foo"), MockGap("bar")]
		writer = MockWriter()

		GapUtility(original_read_len, delimiter).write_gap_file(gaps, writer, additional_header_lines)

		actual_lines = writer.lines()
		self.assertEqual(["#hoopy", "#frood", Gap.header(), "foo", "bar"], actual_lines)

def init_gap(chromosome, gap_start, split_read_name):
	return Gap(split_read_name, chromosome, 0, gap_start, 16, 64)

class MockGap():
	def __init__(self, format_string):
		self._format_string = format_string

	def _format(self, delimiter):
		return self._format_string


class MockWriter():
	def __init__(self):
		self._content = []

	def write(self, content):
		self._content.append(content)
		
	def lines(self):
		return "".join(self._content).splitlines()

if __name__ == "__main__":
	unittest.main() 
