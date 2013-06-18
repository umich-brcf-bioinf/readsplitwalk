#! /usr/bin/env python

import sys ; sys.path.insert(0, "../bin")
import unittest
from IdentifyPairs import BowtieSplitReadBuilder, LegacySplitReadBuilder, SplitRead, _build_read_groups, _write_pairs, _build_pairs_from_groups, _identify_common_group_keys

class LegacySplitReadBuilderTest(unittest.TestCase):

	def test_build(self):
		read_len = 30		
		builder = LegacySplitReadBuilder(read_len, "|")
		split_read = builder.build("name|L|10|strand|chr|100|seq|quality|5")
	
		self.assertEqual("name", split_read.name)
		self.assertEqual("L", split_read.side)
		self.assertEqual(10, split_read.split_len)
		self.assertEqual("strand", split_read._strand)
		self.assertEqual("chr", split_read._chr)
		self.assertEqual(100, split_read._position)
		self.assertEqual(5, split_read._matches)

	def test_build_raisesOnMalformedInput(self):
		builder = LegacySplitReadBuilder(30,"|")
		self.assertRaises(Exception, builder.build, "name|L|10")
		


class BowtieSplitReadBuilderTest(unittest.TestCase):

	def test_build(self):
		read_len = 30		
		builder = BowtieSplitReadBuilder(read_len, "|")
		split_read = builder.build("hw1-name-L-10|strand|chr|100|seq|quality|5")
	
		self.assertEqual("hw1-name", split_read.name)
		self.assertEqual("L", split_read.side)
		self.assertEqual(10, split_read.split_len)
		self.assertEqual("strand", split_read._strand)
		self.assertEqual("chr", split_read._chr)
		self.assertEqual(100, split_read._position)
		self.assertEqual(5, split_read._matches)

	def test_build_raisesOnMalformedInput(self):
		builder = BowtieSplitReadBuilder(30,"|")
		self.assertRaises(Exception, builder.build, "name|L|10")
		
		
class SplitReadTest(unittest.TestCase):
		
	def test_side(self):
		sr = SplitRead(**initParams({'side':"L"}))
		self.assertEqual("L", sr.side)

	def test_name(self):
		sr = SplitRead(**initParams({'name':"the_read_name"}))
		self.assertEqual("the_read_name", sr.name)
		
	def test_format(self):
		params = initParams({'name':"foo", 'side':"L", 'split_len':10, 'strand':"strand", 'chr':"chr", 'position':100, 'matches':5}) 
		sr = SplitRead(**params)
		self.assertEqual("foo~L~10~strand~chr~100~5", sr.format("~"))	

	def test_distance(self):
		srL = SplitRead(**initParams({'position':25}))
		srR = SplitRead(**initParams({'position':100}))
		self.assertEqual(75, srL.distance(srR))
		self.assertEqual(75, srR.distance(srL))
		
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
	

class IdentifyPairsTest(unittest.TestCase):

	def test_identify_common_group_keys(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key1", "R")
		read_len = read1.split_len + read2.split_len
		builder = MockSplitReadBuilder({'read1': read1, 'read2' : read2})
		reader = ["read1", "read2"]

		group_keys = _identify_common_group_keys(builder, reader, MockLogger(), read_len)
	
		self.assertEqual(1, len(group_keys))
		self.assertEqual(True, "key1" in group_keys)

	def test_identify_common_group_keys_noCommonKeys(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key2", "R")
		read_len = read1.split_len + read2.split_len
		builder = MockSplitReadBuilder({'read1': read1, 'read2' : read2})
		reader = ["read1", "read2"]
		group_keys = _identify_common_group_keys(builder, reader, MockLogger(), read_len)
		self.assertEqual(0, len(group_keys))

	def test_identify_common_group_keys_keyOnlyOnOneSide(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key1", "L")
		read_len = read1.split_len + read2.split_len
		builder = MockSplitReadBuilder({'read1': read1, 'read2' : read2})
		reader = ["read1", "read2"]		
		group_keys = _identify_common_group_keys(builder, reader, MockLogger(), read_len)
		self.assertEqual(0, len(group_keys))

	def test_build_read_groups_singleRead(self):
		read1 = MockSplitRead("key1","L")
		builder = MockSplitReadBuilder({'read1' : read1 })
		reader = ["read1"]
		common_keys = set(["key1"])
		
		pairs = _build_read_groups(common_keys, builder, reader, MockLogger())
		
		self.assertEqual(1, len(pairs))
		self.assertEqual(([read1],[]), pairs["key1"])

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

        def test_build_read_groups_filtersOnCommonKey(self):
                read1 = MockSplitRead("key1", "L")
                read2 = MockSplitRead("key2", "R")
                split_read_builder = MockSplitReadBuilder({'read1':read1, 'read2':read2})
                reader = ["read1", "read2"]
                common_keys = set(["key1"])

                pairs = _build_read_groups(common_keys, split_read_builder, reader, MockLogger())

                self.assertEqual(1, len(pairs))
                self.assertEqual(([read1],[]), pairs["key1"])


	def test_build_read_groups_twoLeftReads(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key1", "L")
		split_read_builder = MockSplitReadBuilder({'read1':read1, 'read2':read2})
		reader = ["read1", "read2"]
		common_keys = set(["key1"])	
	
		pairs = _build_read_groups(common_keys, split_read_builder, reader, MockLogger())
		
		self.assertEqual(1, len(pairs))
		self.assertEqual(([read1, read2],[]), pairs["key1"])

	def test_build_read_groups_twoLeftReadsAndTwoRight(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key1", "L")
		read3 = MockSplitRead("key1", "R")
		read4 = MockSplitRead("key1", "R")
		split_read_builder = MockSplitReadBuilder({'1':read1, '2':read2, '3':read3, '4':read4})
		reader = ["1","2","3","4"]
		common_keys = set(["key1"])	

		pairs = _build_read_groups(common_keys, split_read_builder, reader, MockLogger())
		
		self.assertEqual(1, len(pairs))
		self.assertEqual(([read1, read2],[read3, read4]), pairs["key1"])

	def test_write_pairs(self):
		writer = MockWriter()
		leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead")
		rightA = MockSplitRead("key1", "L", "rightA", "rightFormattedRead")		
		pairs=[(leftA, rightA, 5)]
		
		_write_pairs(pairs, writer, MockLogger(), "|")	
		
		self.assertEqual(1, len(writer.lines()))
		self.assertEqual(["leftFormattedRead|rightFormattedRead|5"], writer.lines())
	

	def test_build_pairs_from_groups_simpleDiad(self):
		leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 5)
		rightA = MockSplitRead("key1", "L", "rightA", "rightFormattedRead", 42)		
		groups = {'1':([leftA],[rightA])}
		
		actual = _build_pairs_from_groups(groups, MockLogger())	
		
		expected=[(leftA, rightA, 5)]
		self.assertEqual(1, len(actual))
		self.assertEqual(expected, actual)


	def test_build_pairs_from_groups(self):
		leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 5)
		leftB = MockSplitRead("key1", "L", "leftB", "leftFormattedRead", 10)
		rightA = MockSplitRead("key1", "R", "rightA", "rightFormattedRead", 42)		
		rightB = MockSplitRead("key1", "R", "rightB", "rightFormattedRead", 42)		

		groups = {'1':([leftA, leftB], [rightA, rightB])}
		
		actual = _build_pairs_from_groups(groups, MockLogger())	
	
		actual = sorted(actual)	
		expected=sorted([(leftA, rightA, 5), (leftA, rightB, 5), (leftB, rightA, 10), (leftB, rightB, 10)])
	
		self.assertEqual(4, len(actual))
		self.assertEqual(expected, actual)
		
	
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
		self.name = name
		self.side = side
		self._format = format
		self._distance = distance
		self.split_len = split_len

	def format(self, delimiter="\t"):
		return self._format
		
	def distance(self, distance):
		return self._distance
		
	def key(self):
		return self._key

	def __repr__(self):
		return self.name


class MockSplitReadBuilder():
	def __init__(self, split_reads):
		self._split_reads = split_reads
	
	def build(self, line):
		return self._split_reads[line]



def initParams(updates):
	params = {'name':"name", 'side':"L", 'split_len':10, 'strand':"strand", 'chr':"chr", 'position':100, 'matches':5, 'read_len': 33} 
	params.update(updates);
	return params



if __name__ == "__main__":
    unittest.main() 
