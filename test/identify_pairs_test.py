import unittest
from bin.identify_pairs import BowtieSplitReadBuilder, LegacySplitReadBuilder, ReadLengthValidator, ReadLengthValidationError, SamSplitReadBuilder, SplitRead, _build_read_groups, _write_rsw_pairs, _write_sam_pairs, _build_pairs_from_groups, _identify_common_group_keys, _filter_pairs, _distance_filter, _orientation_filter, _composite_filter


class LegacySplitReadBuilderTestCase(unittest.TestCase):

    def test_build(self):
        read_len = 30       
        builder = LegacySplitReadBuilder(read_len, "|")
        split_read = builder.build("name|L|10|strand|chr|100|seq|quality|5")
    
        self.assertEqual("name", split_read._name)
        self.assertEqual("L", split_read._side)
        self.assertEqual(10, split_read._split_len)
        self.assertEqual("strand", split_read._strand)
        self.assertEqual("chr", split_read._chr)
        self.assertEqual(100, split_read._position)
        self.assertEqual(5, split_read._matches)

    def test_build_raisesOnMalformedInput(self):
        builder = LegacySplitReadBuilder(30,"|")
        self.assertRaises(Exception, builder.build, "name|L|10")
        


class BowtieSplitReadBuilderTestCase(unittest.TestCase):

    def test_build(self):
        read_len = 30       
        builder = BowtieSplitReadBuilder(read_len, "|")
        split_read = builder.build("hw1-name-L-10|strand|chr|100|seq|quality|5")
    
        self.assertEqual("hw1-name", split_read._name)
        self.assertEqual("L", split_read._side)
        self.assertEqual(10, split_read._split_len)
        self.assertEqual("strand", split_read._strand)
        self.assertEqual("chr", split_read._chr)
        self.assertEqual(100, split_read._position)
        self.assertEqual(5, split_read._matches)

    def test_build_raisesOnMalformedInput(self):
        builder = BowtieSplitReadBuilder(30,"|")
        self.assertRaises(Exception, builder.build, "name|L|10")
        
        
class SamSplitReadBuilderTestCase(unittest.TestCase):

    def test_build(self):
        read_len = 30       
        builder = SamSplitReadBuilder(read_len, "|")
        split_read = builder.build("hw1:name-L-10|16|chr|100|25|cigar|*|0|0|GCAGT|DDDCC@|NM:i:0  X0:i:1")
    
        self.assertEqual("hw1:name", split_read._name)
        self.assertEqual("L", split_read._side)
        self.assertEqual(10, split_read._split_len)
        self.assertEqual("-", split_read._strand)
        self.assertEqual("chr", split_read._chr)
        self.assertEqual(100, split_read._position)
        self.assertEqual(None, split_read._matches)

    def test_build_raisesOnMalformedInput(self):
        builder = BowtieSplitReadBuilder(30,"|")
        self.assertRaises(Exception, builder.build, "name|L|10")
        
        
class SplitReadTestCase(unittest.TestCase):
        
    def test_side(self):
        sr = SplitRead(**initParams({'side':"L"}))
        self.assertEqual("L", sr._side)

    def test_name(self):
        sr = SplitRead(**initParams({'name':"the_read_name"}))
        self.assertEqual("the_read_name", sr._name)
        
    def test_format(self):
        params = initParams({'name':"foo", 'side':"L", 'split_len':10, 'strand':"strand", 'chromosome':"chr", 'position':100, 'matches':5}) 
        sr = SplitRead(**params)
        self.assertEqual("foo~L~10~strand~chr~100~5", sr.format("~"))   

    def test_gap_distance(self):
        srL1 = SplitRead(**initParams({'position':25, 'split_len':60, 'strand':"+", 'side':"L", 'original_read_len':200}))
        srR1 = SplitRead(**initParams({'position':100, 'split_len':40, 'strand':"+", 'side':"R", 'original_read_len':200}))
        self.assertEqual(15, srL1.gap_distance(srR1))
        self.assertEqual(15, srR1.gap_distance(srL1))
        
        srL2 = SplitRead(**initParams({'position':100, 'split_len':40, 'strand':"-", 'side':"L", 'original_read_len':200}))
        srR2 = SplitRead(**initParams({'position':55, 'split_len':30, 'strand':"-", 'side':"R", 'original_read_len':200}))
        self.assertEqual(15, srL2.gap_distance(srR2))
        self.assertEqual(15, srR2.gap_distance(srL2))
        
    def test_key_leftKeyPassesThrough(self):
        read_len = 30
        builder = LegacySplitReadBuilder(read_len, "-")
        split_read = builder.build("name-L-10-strand-chr-100-seq-quality-5-foo-bar\n")
        self.assertEqual("name|L|10|strand|chr", split_read.key())

    def test_key_rightKeySwitchesSideAndSplitLength(self):
        read_len = 30
        builder = LegacySplitReadBuilder(read_len, "-")
        split_read = builder.build("name-R-20-strand-chr-100-seq-quality-5-foo-bar\n")
        self.assertEqual("name|L|10|strand|chr", split_read.key())
    
    def test_is_oriented_true(self):
        left = SplitRead(**initParams({'side':"L", 'position': 1, 'strand':"+"}))
        right = SplitRead(**initParams({'side':"R", 'position': 2, 'strand':"+"}))
        self.assertEqual(True, left.is_oriented(right))
        self.assertEqual(True, right.is_oriented(left))

        left = SplitRead(**initParams({'side':"L", 'position': 2, 'strand':"-"}))
        right = SplitRead(**initParams({'side':"R", 'position': 1, 'strand':"-"}))
        self.assertEqual(True, left.is_oriented(right))
        self.assertEqual(True, right.is_oriented(left))

    def test_is_oriented_false(self):
        left = SplitRead(**initParams({'side':"L", 'position': 2, 'strand':"+"}))
        right = SplitRead(**initParams({'side':"R", 'position': 1, 'strand':"+"}))
        self.assertEqual(False, left.is_oriented(right))
        self.assertEqual(False, right.is_oriented(left))

        left = SplitRead(**initParams({'side':"L", 'position': 1, 'strand':"-"}))
        right = SplitRead(**initParams({'side':"R", 'position': 2, 'strand':"-"}))
        self.assertEqual(False, left.is_oriented(right))
        self.assertEqual(False, right.is_oriented(left))

    def test_is_oriented_falseIfEqualPositions(self):
        left = SplitRead(**initParams({'side':"L", 'position': 1, 'strand':"+"}))
        right = SplitRead(**initParams({'side':"R", 'position': 1, 'strand':"+"}))
        self.assertEqual(False, left.is_oriented(right))
        self.assertEqual(False, right.is_oriented(left))

    def test_is_oriented_falseIfStrandDifferent(self):
        left = SplitRead(**initParams({'side':"L", 'position': 10, 'strand':"+"}))
        right = SplitRead(**initParams({'side':"R", 'position': 150, 'strand':"-"}))
        self.assertEqual(False, left.is_oriented(right))
        self.assertEqual(False, right.is_oriented(left))

    def test_is_oriented_falseIfSideSame(self):
        left = SplitRead(**initParams({'side':"L", 'position': 10, 'strand':"+"}))
        left2 = SplitRead(**initParams({'side':"L", 'position': 150, 'strand':"+"}))
        self.assertEqual(False, left.is_oriented(left2))
        self.assertEqual(False, left2.is_oriented(left))

    def test_left_name(self):
        left = SplitRead(**initParams({'name':'readA', 'side':"L", 'split_len': 10, 'original_read_len': 100}))
        right = SplitRead(**initParams({'name':'readA', 'side':"R", 'split_len': 90, 'original_read_len': 100}))
        self.assertEqual("readA-L-10", left.left_name())
        self.assertEqual("readA-L-10", right.left_name())
        

    def test_write_sam_pairs_skipsOrphanedLines(self):
        writer = MockWriter()
        
        leftA = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 10}))
        rightA = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 5}))
        
        read_group_key = leftA.key()
        read_group_pairs = {read_group_key:[(leftA, rightA, 5)]}

        split_read_from_file = SplitRead(**initParams({'name':'readB', 'side':"L", 'position': 10}))
        split_read_from_file.write_sam_pairs(read_group_pairs, "line", writer)
        
        self.assertEqual(0, len(writer.lines()))

    def test_write_sam_pairs_skipsIfWrongPosition(self):
        writer = MockWriter()
        
        leftA = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 10}))
        rightA = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 5}))
        
        read_group_key = leftA.key()
        read_group_pairs = {read_group_key:[(leftA, rightA, 5)]}

        split_read_from_file = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 999}))
        split_read_from_file.write_sam_pairs(read_group_pairs, "line\t"*12, writer)
        
        self.assertEqual(0, len(writer.lines()))
        
    def test_write_sam_pairs_writesLineForEachPairParticipation(self):
        writer = MockWriter()
        stub_line = "readA|"*12
        leftA5 = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 5}))
        rightA10 = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 10}))
        rightA15 = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 15}))
        
        read_group_key = leftA5.key()
        read_group_pairs = {read_group_key:[(leftA5, rightA10), (leftA5, rightA15)]}

        split_read_from_file = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 5}))
        split_read_from_file.write_sam_pairs(read_group_pairs, stub_line+"\n", writer, "|")
        
        actual_lines = writer.lines()
        self.assertEqual(2, len(actual_lines))

    def test_write_sam_pairs_writesPositiveStrandSamLines(self):
        writer = MockWriter()
        input_line = "readA 147 chr12   5   255 42M *   0   0   TCACC   DDDDD   XA:i:0"
        left_read = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 5000, 'split_len':15, 'original_read_len':100}))
        right_read = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 5200, 'split_len':85, 'original_read_len':100}))
        
        read_group_key = left_read.key()
        read_group_pairs = {read_group_key:[(left_read, right_read)]}

        first_split_read_from_file = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 5000, 'split_len':15, 'original_read_len':100}))
        first_split_read_from_file.write_sam_pairs(read_group_pairs, input_line+"\n", writer)
        second_split_read_from_file = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 5200, 'split_len':85, 'original_read_len':100}))
        second_split_read_from_file.write_sam_pairs(read_group_pairs, input_line+"\n", writer)
        
        actual_lines = writer.lines()
        self.assertEqual(2, len(actual_lines))
        self.assertEqual("readA-L-15    67  chr12   5000    255 42M =   5200    200 TCACC   DDDDD   XA:i:0", actual_lines[0])
        self.assertEqual("readA-L-15    131 chr12   5200    255 42M =   5000    -200    TCACC   DDDDD   XA:i:0", actual_lines[1])

    def test_write_sam_pairs_writesNegativeStrandSamLines(self):
        writer = MockWriter()
        input_line = "readA 147 chr12   5   255 42M *   0   0   TCACC   DDDDD   XA:i:0"
        left_read = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 5200, 'split_len':15, 'original_read_len':100}))
        right_read = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 5000, 'split_len':85, 'original_read_len':100}))
        
        read_group_key = left_read.key()
        read_group_pairs = {read_group_key:[(left_read, right_read)]}

        first_split_read_from_file = SplitRead(**initParams({'name':'readA', 'side':"L", 'position': 5200, 'split_len':15, 'original_read_len':100}))
        first_split_read_from_file.write_sam_pairs(read_group_pairs, input_line+"\n", writer)
        second_split_read_from_file = SplitRead(**initParams({'name':'readA', 'side':"R", 'position': 5000, 'split_len':85, 'original_read_len':100}))
        second_split_read_from_file.write_sam_pairs(read_group_pairs, input_line+"\n", writer)
        
        actual_lines = writer.lines()
        self.assertEqual(2, len(actual_lines))
        self.assertEqual("readA-L-15    115 chr12   5200    255 42M =   5000    -200    TCACC   DDDDD   XA:i:0", actual_lines[0])
        self.assertEqual("readA-L-15    147 chr12   5000    255 42M =   5200    200 TCACC   DDDDD   XA:i:0", actual_lines[1])



    def test_add_to_read_groups_doesNothingWhenNotInCommonKeys(self):
        readA = SplitRead(**initParams({'name':'readA'}))
        common_keys = set([readA.key()])
        read_groups = {}

        readB = SplitRead(**initParams({'name':'readB'}))       
        readB.add_to_read_groups(common_keys, read_groups)
        
        self.assertEqual(0, len(read_groups))

    def test_add_to_read_groups_addsToCorrectSide(self):
        left10 = SplitRead(**initParams({'name':'readA', 'side':"L", 'position':10, 'split_len':40, 'original_read_len':100}))
        left15 = SplitRead(**initParams({'name':'readA', 'side':"L", 'position':15, 'split_len':40, 'original_read_len':100}))
        right30 = SplitRead(**initParams({'name':'readA', 'side':"R", 'position':30, 'split_len':60, 'original_read_len':100}))
        read_group_key = left10.key()
        common_keys = set([read_group_key])
        read_groups = {}

        left10.add_to_read_groups(common_keys, read_groups)
        left15.add_to_read_groups(common_keys, read_groups)
        right30.add_to_read_groups(common_keys, read_groups)
        
        self.assertEqual(1, len(read_groups))
        self.assertEqual(([left10, left15], [right30]), read_groups[read_group_key])

    def test_add_to_group_keys(self):
        group_keys = {'L':set(), 'R':set()}
        readA10 = SplitRead(**initParams({'name':'readA', 'side':"L", 'position':10}))
        readA15 = SplitRead(**initParams({'name':'readA', 'side':"L", 'position':15}))
        readB30 = SplitRead(**initParams({'name':'readB', 'side':"R", 'position':30}))
        
        readA10.add_to_group_keys(group_keys)
        readA15.add_to_group_keys(group_keys)
        readB30.add_to_group_keys(group_keys)
        
        self.assertEqual(2, len(group_keys))
        self.assertEqual(set([readA10.key(), readA15.key()]), group_keys["L"])
        self.assertEqual(set([readB30.key()]), group_keys["R"])


class ReadLengthValidatorTestCase(unittest.TestCase):
    
    def test_check_split_length(self):
        read_length = 42
        validator = ReadLengthValidator(read_length)
        validator.check_split_length(42)
        self.assertRaises(ReadLengthValidationError, validator.check_split_length, 43)

    def test_check_read_length(self):
        read_length = 50
        validator = ReadLengthValidator(read_length)
        validator.check_split_length(40)
        validator.check_split_length(10)
        validator.check_read_length()

        validator = ReadLengthValidator(read_length)
        validator.check_split_length(40)
        validator.check_split_length(11)
        self.assertRaises(ReadLengthValidationError, validator.check_read_length)

        validator = ReadLengthValidator(read_length)
        validator.check_split_length(15)
        validator.check_split_length(10)
        self.assertRaises(ReadLengthValidationError, validator.check_read_length)


class IdentifyPairsTestCase(unittest.TestCase):

    def test_identify_common_group_keys(self):
        read1 = MockSplitRead("key1", "L")
        read2 = MockSplitRead("key1", "R")
        builder = MockSplitReadBuilder({'read1': read1, 'read2': read2}, ["@header1"])
        reader = ["@header1", "read1", "read2"]

        group_keys = _identify_common_group_keys(builder, MockValidator(), reader, MockLogger())
    
        self.assertEqual(1, len(group_keys))
        self.assertEqual(True, "key1" in group_keys)

    def test_identify_common_group_keys_noCommonKeys(self):
        read1 = MockSplitRead("key1", "L")
        read2 = MockSplitRead("key2", "R")
        read_len = read1.split_len + read2.split_len
        builder = MockSplitReadBuilder({'read1': read1, 'read2' : read2})
        reader = ["read1", "read2"]

        group_keys = _identify_common_group_keys(builder, MockValidator(), reader, MockLogger())

        self.assertEqual(0, len(group_keys))

    def test_identify_common_group_keys_keyOnlyOnOneSide(self):
        read1 = MockSplitRead("key1", "L")
        read2 = MockSplitRead("key1", "L")
        read_len = read1.split_len + read2.split_len
        builder = MockSplitReadBuilder({'read1': read1, 'read2' : read2})
        reader = ["read1", "read2"]     

        group_keys = _identify_common_group_keys(builder, MockValidator(), reader, MockLogger())

        self.assertEqual(0, len(group_keys))

    def test_build_read_groups_twoDistinctReads(self):
        read1 = MockSplitRead("key1", "L")
        read2 = MockSplitRead("key2", "R")
        split_read_builder = MockSplitReadBuilder({'read1':read1, 'read2':read2})
        reader = ["read1", "read2"]
        common_keys = set(["key1","key2"])  
    
        pairs = _build_read_groups(common_keys, split_read_builder, reader, MockLogger())
        
        self.assertEqual(2, len(pairs))
        self.assertEqual(([read1],[]), pairs["key1"])
        self.assertEqual(([],[read2]), pairs["key2"])

    def test_build_read_groups_skipsHeaderLines(self):
        read1 = MockSplitRead("key1", "L")
        read2 = MockSplitRead("key1", "L")
        split_read_builder = MockSplitReadBuilder({'1':read1, '2':read2}, ["@h1", "@h2"])
        reader = ["@h1", "@h2", "1","2"]
        common_keys = set() 

        pairs = _build_read_groups(common_keys, split_read_builder, reader, MockLogger())
        
        self.assertEqual(1, len(pairs))

    def test_write_rsw_pairs(self):
        writer = MockWriter()
        leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 5)
        rightA = MockSplitRead("key1", "L", "rightA", "rightFormattedRead")     
        pairs = {'key1': [(leftA, rightA)]}
        
        _write_rsw_pairs(pairs, writer, MockLogger(), "|")  
        
        self.assertEqual(1, len(writer.lines()))
        self.assertEqual(["leftFormattedRead|rightFormattedRead|5"], writer.lines())


    def test_write_sam_pairs_headersPassThrough(self):
        reader = ["@header1\n","@header2\n"]
        builder = MockSplitReadBuilder({}, set(["@header1", "@header2"]))
        writer = MockWriter()
        read_group_pairs = {}
        
        _write_sam_pairs(read_group_pairs, reader, builder, writer, MockLogger())   
        
        actual_lines = writer.lines()
        self.assertEqual(2, len(actual_lines))
        self.assertEqual("@header1", actual_lines[0])
        self.assertEqual("@header2", actual_lines[1])

    def test_write_sam_pairs(self):
        reader = ["line1\n", "line2\n"]
        read1 = MockSplitRead("key1", "L")
        read2 = MockSplitRead("key2", "R")
        builder = MockSplitReadBuilder({'line1': read1, 'line2' : read2}, set(["@header1", "@header2"]))
        read_group_pairs = {}
        
        _write_sam_pairs(read_group_pairs, reader, builder, MockWriter(), MockLogger()) 
        
        self.assertEqual(1, read1.write_sam_pairs_called)
        self.assertEqual(1, read2.write_sam_pairs_called)

    def test_build_pairs_from_groups_simpleDiad(self):
        leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead")
        rightA = MockSplitRead("key1", "L", "rightA", "rightFormattedRead")     
        groups = {'key1':([leftA],[rightA])}
        
        actual_pairs = _build_pairs_from_groups(groups, MockLogger()) 
        expected=[(leftA, rightA)]
        self.assertEqual(1, len(actual_pairs))
        self.assertEqual(expected, actual_pairs['key1'])

    def test_build_pairs_from_groups(self):
        leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead")
        leftB = MockSplitRead("key1", "L", "leftB", "leftFormattedRead")
        rightA = MockSplitRead("key1", "R", "rightA", "rightFormattedRead")     
        rightB = MockSplitRead("key1", "R", "rightB", "rightFormattedRead")     

        groups = {'key1':([leftA, leftB], [rightA, rightB])}
        
        actual_pairs = _build_pairs_from_groups(groups, MockLogger()) 
    
        self.assertEqual(1, len(actual_pairs))
        sorted_actual_pairs = sorted(actual_pairs['key1']) 
        sorted_expected_pairs = sorted([(leftA, rightA), (leftA, rightB), (leftB, rightA), (leftB, rightB)])
    
        self.assertEqual(4, len(sorted_actual_pairs))
        self.assertEqual(sorted_expected_pairs, sorted_actual_pairs)
        

    def test_filter_pairs(self):
        input_pairs = {'key1':[("READ1", "read2"), ("read3","read4")], 'key2':[("READ5", "read6")]}
        def mock_filter(read1, read2): return read1 == read1.upper()

        actual_pairs = _filter_pairs(input_pairs, mock_filter, MockLogger())

        self.assertEqual(2, len(actual_pairs))
        self.assertEqual([("READ1", "read2")], actual_pairs["key1"])
        self.assertEqual([("READ5", "read6")], actual_pairs["key2"])


    def test_distance_filter(self):
        filter = _distance_filter(5, 10)
        leftA4 = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 4)
        leftA5 = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 5)
        leftA7 = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 7)
        leftA10 = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 10)
        leftA11 = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 11)
        rightA = MockSplitRead("key1", "R", "rightA", "rightFormattedRead", 42)

        self.assertEqual(False, filter(leftA4, rightA))
        self.assertEqual(True, filter(leftA5, rightA))
        self.assertEqual(True, filter(leftA7, rightA))
        self.assertEqual(True, filter(leftA10, rightA))
        self.assertEqual(False, filter(leftA11, rightA))

        
    def test_orientation_filter(self):
        left1 = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 4)
        left1._is_oriented = False
        left2 = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 5)
        left2._is_oriented = True
        right = MockSplitRead("key1", "R", "rightA", "rightFormattedRead", 42)

        self.assertEqual(False, _orientation_filter(left1, right))
        self.assertEqual(True, _orientation_filter(left2, right))

    def test_composite_filter(self):
        def filter1(read1, read2):
            return read1 == read1.upper()
        def filter2(read1, read2):
            return read2 == read2.lower()
        composite_filter = _composite_filter([filter1, filter2])
        
        self.assertEqual(True, composite_filter("READ", "read"))            
        self.assertEqual(False, composite_filter("Read", "read"))           
        self.assertEqual(False, composite_filter("READ", "Read"))           
        
        
class MockFilter():

    def __init__(self):
        MockFilter.filter_pair_was_called = 0

    def filter_pair(self, read1, read2):
        MockFilter.filter_pair_was_called += 1

    
class MockLogger():
    def log(self, message):
        pass
    
    def log_usage(self, message):
        pass        

class MockWriter():
    def __init__(self):
        self._content = []

    def write(self, content):
        self._content.append(content)
        
    def lines(self):
        return "".join(self._content).splitlines()


class MockPairRowFormatter():
    def format(self, leftRead, rightRead):
        return "rowHeader"+leftRead[-1]+rightRead[-1]+"|"+leftRead+"|"+rightRead+"|rowFooter"+leftRead[-1]+rightRead[-1]+"\n"

        
class MockSplitRead():
    def __init__(self, key, side, name = "name", format = "format", distance=42, split_len=33):
        
        self._key = key
        self._name = name
        self._side = side
        self._format = format
        self._distance = distance
        self.split_len = split_len
        self.write_sam_pairs_called = 0
        self._is_oriented = False

    def format(self, delimiter="\t"):
        return self._format
        
    def gap_distance(self, distance):
        return self._distance
        
    def key(self):
        return self._key

    def __repr__(self):
        return self._name

    def is_oriented(self, other):       
        return self._is_oriented

    def write_sam_pairs(self, read_group_pairs, line, writer, delimiter):
        self.write_sam_pairs_called += 1

    def add_to_read_groups(self, common_keys, read_groups):
        read_group = read_groups.setdefault(self._key, ([], [])) 
        index = 0 if self._side == 'L' else 1
        read_group[index].append(self)
        
    def add_to_group_keys(self, group_keys):
        group_keys[self._side].add(self._key)
        
    def check_split_length(self, validator): 
        pass

class MockSplitReadBuilder():
    def __init__(self, split_reads, headers=set()):
        self._split_reads = split_reads
        self._headers = headers 
    def build(self, line):
        return self._split_reads[line.rstrip()]

    def is_header(self, line):
        return line.rstrip() in self._headers

class MockValidator():
    def check_split_length(self, split_len): pass

    def check_read_length(self): pass

def initParams(updates):
    params = {'name':"name", 'side':"L", 'split_len':10, 'strand':"+", 'chromosome':"chr", 'position':100, 'matches':5, 'original_read_len': 33} 
    params.update(updates);
    return params



if __name__ == "__main__":
    unittest.main() 
