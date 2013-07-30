import unittest
from bin.cluster_gaps import GapUtility, Gap, MissingReadGroupError, InvalidReadGroupError


class GapTestCase(unittest.TestCase):

    def test_gap_width(self):
        self.assertEqual(12, Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, 16, 64).gap_width())

    def test_read_width(self):
        self.assertEqual(64, Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, 16, 64)._read_width())

    def test_original_read_name(self):
        self.assertEqual("split-read-name", Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, 16, 64)._original_read_name())

    def test_header(self):
        self.assertEqual("#chromosome|cluster|sample_name|gap_start|gap_end|gap_width|read_start|read_end|read_width|split_read_name|original_read_name", Gap.header("|"))
        
    def test_format(self):
        gap = Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, 16, 64)
        self.assertEqual("chromosome|-1|sampleName|4|16|12|0|64|64|split-read-name-L-13|split-read-name", gap.format("|"))
        gap.cluster = "42"
        self.assertEqual("chromosome|42|sampleName|4|16|12|0|64|64|split-read-name-L-13|split-read-name", gap.format("|"))

    def test_eq(self):
        base = Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, 16, 64)
        self.assertEqual(True, base == Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, 16, 64))
        self.assertEqual(False, base == "foo")
        self.assertEqual(False, base == Gap("X", "split-read-name-L-13", "chromosome", 0, 4, 16, 64))
        self.assertEqual(False, base == Gap("sampleName", "X", "chromosome", 0, 4, 16, 64))
        self.assertEqual(False, base == Gap("sampleName", "split-read-name-L-13", "chromosome", -42, 4, 16, 64))
        self.assertEqual(False, base == Gap("sampleName", "split-read-name-L-13", "chromosome", 0, -42, 16, 64))
        self.assertEqual(False, base == Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, -42, 64))
        self.assertEqual(False, base == Gap("sampleName", "split-read-name-L-13", "chromosome", 0, 4, 16, -42))


class GapUtilityTestCase(unittest.TestCase):

    def test_build_gap_leftmost(self):   
        gap_utility = GapUtility(original_read_len=50, delimiter="|", logger=MockLogger())
        gap_utility._read_group_sample_dict = {'1':'sampleName'}
        sam_line = "read-name-L-42|67|transcript42|100|score|cigar|=|150|50|ACGCT|CCCFF|RG:Z:1"

        gap = gap_utility.build_gap(sam_line)
        
        self.assertEqual("transcript42", gap.chromosome)
        self.assertEqual(105, gap.gap_start)
        self.assertEqual(150, gap._gap_end)
        self.assertEqual(100, gap._read_start)
        self.assertEqual(195, gap._read_end)
        self.assertEqual("read-name-L-42", gap._split_read_name)

    def test_build_gap_rightmost(self):   
        gap_utility = GapUtility(original_read_len=50, delimiter="|", logger=MockLogger())
        gap_utility._read_group_sample_dict = {'1':'sampleName'}
        sam_line = "read-name-L-42|131|transcript42|150|score|cigar|=|80|-50|ACGCT|CCCFF|RG:Z:1"

        gap = gap_utility.build_gap(sam_line)
        
        self.assertEqual("transcript42", gap.chromosome)
        self.assertEqual(80, gap._read_start)
        self.assertEqual(125, gap.gap_start)
        self.assertEqual(150, gap._gap_end)
        self.assertEqual(155, gap._read_end)
        self.assertEqual("read-name-L-42", gap._split_read_name)


    def test_samfile_to_gaps_ignoresRightmostAlignment(self):     
        original_read_len = 50
        delimiter="|"
        sam_file = ["@header1", "@RG|ID:1|SM:sampleName", 
                "read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|RG:Z:1", 
                "read1-L-1|131|transcript42|200|score|cigar|=|150|-50|ACGCT|qual|RG:Z:1",
                "read2-L-1|131|transcript43|205|score|cigar|=|155|-50|ACGCT|qual|RG:Z:1",
                "read2-L-1|67|transcript43|155|score|cigar|=|205|50|ACGCT|qual|RG:Z:1"]

        gaps = GapUtility(original_read_len, delimiter, MockLogger()).samfile_to_gaps(sam_file)

        self.assertEqual(2, len(gaps))
        self.assertEqual("transcript42", gaps[0].chromosome)
        self.assertEqual("transcript43", gaps[1].chromosome)

    def test_samfile_to_gaps_capturesReadgroupSampleMap(self):     
        original_read_len = 50
        delimiter="|"
        sam_file = ["@RG|ID:foo|SM:adam", "@RG|ID:bar|SM:betty", 
                "read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|stuff1:i:1|RG:Z:foo|stuff2:i:1", 
                "read2-L-1|67|transcript43|155|score|cigar|=|205|50|ACGCT|qual|stuff:i:1|stuff2:i:1|RG:Z:bar"
                ]

        gaps = GapUtility(original_read_len, delimiter, MockLogger()).samfile_to_gaps(sam_file)

        self.assertEqual(2, len(gaps))
        self.assertEqual("adam", gaps[0].sample)
        self.assertEqual("betty", gaps[1].sample)

    def test_samfile_to_gaps_throwsOnInvalidReadGroup(self):     
        original_read_len = 50
        delimiter="|"
        sam_file = ["@RG|ID:foo|SM:adam",
                "read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|stuff1:i:1|RG:Z:argh!|stuff2:i:1"
                ]
        gap_utility = GapUtility(original_read_len, delimiter, MockLogger())
        self.assertRaises(InvalidReadGroupError, gap_utility.samfile_to_gaps, sam_file)

    def test_samfile_to_gaps_throwsOnMissingReadGroup(self):     
        original_read_len = 50
        delimiter="|"
        sam_file = ["@RG|ID:foo|SM:adam",
                "read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual"
                ]
        gap_utility = GapUtility(original_read_len, delimiter, MockLogger())
        self.assertRaises(MissingReadGroupError, gap_utility.samfile_to_gaps, sam_file)

        sam_file = ["@RG|ID:foo|SM:adam",
                "read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|stuff1:i:1|stuff2:i:1"
                ]
        gap_utility = GapUtility(original_read_len, delimiter, MockLogger())
        self.assertRaises(MissingReadGroupError, gap_utility.samfile_to_gaps, sam_file)


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

        actual_gaps = GapUtility(original_read_len, delimiter, MockLogger()).sort_gaps(input_gaps)

        self.assertEqual(len(input_gaps), len(actual_gaps))
        self.assertEqual(sorted_gaps, actual_gaps)

    def test_write_gap_file(self):
        original_read_len = 50
        delimiter = "|"
        additional_header_lines=["hoopy", "frood"]
        gaps = [MockGap("foo"), MockGap("bar")]
        writer = MockWriter()

        GapUtility(original_read_len, delimiter, MockLogger()).write_gap_file(gaps, writer, additional_header_lines)

        actual_lines = writer.lines()
        self.assertEqual(["#hoopy", "#frood", Gap.header("|"), "foo", "bar"], actual_lines)

    def test_write_sam_file(self):
        
        gap_utility = GapUtility(original_read_len=10, delimiter="|", logger=MockLogger())
        gap_utility._read_group_sample_dict = {'1':'sampleName'}

        read1_leftmost = "read1-L-1|67|transcript42|150|score|cigar|=|200|50|ACGCT|qual|RG:Z:1"
        read1_rightmost = "read1-L-1|131|transcript42|200|score|cigar|=|150|-50|GCAGG|qual|RG:Z:1"
        read2_leftmost =  "read2-L-1|147|transcript43|155|score|cigar|=|205|50|ACGCT|qual|RG:Z:1"
        read2_rightmost = "read2-L-1|115|transcript43|205|score|cigar|=|155|-50|GCAGG|qual|RG:Z:1"

        input_sam_file = \
            [line +"\n" for line in ["@header1", "@header2", 
                read1_leftmost, read1_rightmost,
                read2_leftmost, read2_rightmost
                ]
            ]
                    
        read1_gap = gap_utility.build_gap(read1_leftmost)        
        read1_gap.cluster = 5
        read2_gap = gap_utility.build_gap(read2_leftmost)
        read2_gap.cluster = 10
        gaps = [read1_gap, read2_gap]
    
        writer = MockWriter()
        additional_header_lines=["hoopy", "frood"]

        gap_utility.write_sam_file(input_sam_file, gaps, writer, additional_header_lines)
        
        actual_lines = writer.lines()        
        self.assertEqual(["@CO\thoopy", "@CO\tfrood"], actual_lines[0:2])
        self.assertEqual(["@header1", "@header2"], actual_lines[2:4])
        self.assertEqual([read1_leftmost + "|XC:i:5", read1_rightmost + "|XC:i:5"], actual_lines[4:6])
        self.assertEqual([read2_leftmost + "|XC:i:10", read2_rightmost + "|XC:i:10"], actual_lines[6:8])


def init_gap(chromosome, gap_start, split_read_name):
    return Gap("sampleName", split_read_name, chromosome, 0, gap_start, 16, 64)

class MockGap():
    def __init__(self, format_string):
        self._format_string = format_string

    def format(self, delimiter):
        return self._format_string


class MockWriter():
    def __init__(self):
        self._content = []

    def write(self, content):
        self._content.append(content)
        
    def lines(self):
        return "".join(self._content).splitlines()

class MockLogger():
    def log(self, message):
        pass

if __name__ == "__main__":
    unittest.main() 
