#! /usr/bin/env python

import unittest
from IdentifyPairs import SplitReadBuilder, SplitRead, build_read_groups, filter_orphans, print_pairs, build_pairs_from_groups

class SplitReadBuilderTest(unittest.TestCase):

	def test_init(self):
		read_len = 30		
		builder = SplitReadBuilder(read_len, "|")
		split_read = builder.build("name|L|10|strand|chr|100|seq|quality|5")
	
		self.assertEqual(30, split_read._read_len)
		self.assertEqual("name", split_read.name)
		self.assertEqual("L", split_read.side)
		self.assertEqual(10, split_read._split_len)
		self.assertEqual("strand", split_read._strand)
		self.assertEqual("chr", split_read._chr)
		self.assertEqual(100, split_read._position)
		self.assertEqual(5, split_read._matches)

	def test_build_toleratesMissingMatches(self):
		builder = SplitReadBuilder(30,"|")
		self.assertRaises(Exception, builder.build, "name|L|10")


	def test_build_raisesOnMalformedInput(self):
		builder = SplitReadBuilder(30,"|")
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

	def test_key_equalForIdenticalRead(self):
		params = initParams({})
		srA = SplitRead(**params)
		srB = SplitRead(**params)
		self.assertEqual(srA.key(), srB.key())

	def test_key_notEqualForDifferentNameSideLengthStrandChr(self):
		srA = SplitRead(**initParams({'name':"A"}))
		srB = SplitRead(**initParams({'name':"B"}))
		self.assertNotEqual(srA.key(), srB.key())

		srA = SplitRead(**initParams({'side':"L"}))
		srB = SplitRead(**initParams({'side':"R"}))
		self.assertNotEqual(srA.key(), srB.key())

		srA = SplitRead(**initParams({'split_len':10}))
		srB = SplitRead(**initParams({'split_len':20}))
		self.assertNotEqual(srA.key(), srB.key())

		srA = SplitRead(**initParams({'strand':"+"}))
		srB = SplitRead(**initParams({'strand':"-"}))
		self.assertNotEqual(srA.key(), srB.key())

		srA = SplitRead(**initParams({'chr':"chr10"}))
		srB = SplitRead(**initParams({'chr':"chr20"}))
		self.assertNotEqual(srA.key(), srB.key())

	def test_key_equalForDifferentPosition(self):
		srA = SplitRead(**initParams({'position':100}))
		srB = SplitRead(**initParams({'position':110}))
		self.assertEqual(srA.key(), srB.key())

	def test_key_equalsForDifferentMatch(self):
		srA = SplitRead(**initParams({'matches':5}))
		srB = SplitRead(**initParams({'matches':6}))
		self.assertEqual(srA.key(), srB.key())

	def test_key_equalForComplimentaryRead(self):
		srA = SplitRead(**initParams({'read_len':30, 'side':'L', 'split_len':12}))
		srB = SplitRead(**initParams({'read_len':30, 'side':'R', 'split_len':18}))
		self.assertEqual(srA.key(), srB.key())

	def test_distance(self):
		srL = SplitRead(**initParams({'position':25}))
		srR = SplitRead(**initParams({'position':100}))
		self.assertEqual(75, srL.distance(srR))
		self.assertEqual(75, srR.distance(srL))
	

class SplitWalkReadTest(unittest.TestCase):	
	def test_build_read_groups_singleRead(self):
		read1 = MockSplitRead("key1", "L")
		split_read_builder = MockSplitReadBuilder({'read1':read1})
		reader = ["read1"]
		
		pairs = build_read_groups(split_read_builder, reader, MockLogger())
		
		self.assertEqual(1, len(pairs))
		self.assertEqual(([read1],[]), pairs["key1"])

	def test_build_read_groups_twoDistinctReads(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key2", "R")
		split_read_builder = MockSplitReadBuilder({'read1':read1, 'read2':read2})
		reader = ["read1", "read2"]
		
		pairs = build_read_groups(split_read_builder, reader, MockLogger())
		
		self.assertEqual(2, len(pairs))
		self.assertEqual(([read1],[]), pairs["key1"])
		self.assertEqual(([],[read2]), pairs["key2"])

	def test_build_read_groups_twoLeftReads(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key1", "L")
		split_read_builder = MockSplitReadBuilder({'read1':read1, 'read2':read2})
		reader = ["read1", "read2"]
		
		pairs = build_read_groups(split_read_builder, reader, MockLogger())
		
		self.assertEqual(1, len(pairs))
		self.assertEqual(([read1, read2],[]), pairs["key1"])

	def test_build_read_groups_twoLeftReadsAndTwoRight(self):
		read1 = MockSplitRead("key1", "L")
		read2 = MockSplitRead("key1", "L")
		read3 = MockSplitRead("key1", "R")
		read4 = MockSplitRead("key1", "R")
		split_read_builder = MockSplitReadBuilder({'1':read1, '2':read2, '3':read3, '4':read4})
		reader = ["1","2","3","4"]
		
		pairs = build_read_groups(split_read_builder, reader, MockLogger())
		
		self.assertEqual(1, len(pairs))
		self.assertEqual(([read1, read2],[read3, read4]), pairs["key1"])

	def test_filter_orphans(self):
		read_groups={
			'key1':(["left1"], []),
			'key2':(["left2"], ["right2"]),
			'key3':([], ["right3"])
			}
		
		filter_orphans(read_groups, MockLogger())
		
		self.assertEqual(1, len(read_groups))
		self.assertEqual((["left2"], ["right2"]), read_groups["key2"])
		
	def test_print_pairs(self):
		writer = MockWriter()
		leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead")
		rightA = MockSplitRead("key1", "L", "rightA", "rightFormattedRead")		
		pairs=[(leftA, rightA, 5)]
		
		print_pairs(pairs, writer, "|", MockLogger())	
		
		self.assertEqual(1, len(writer.lines()))
		self.assertEqual(["leftFormattedRead|rightFormattedRead|5"], writer.lines())
	

	def test_build_pairs_from_groups_simpleDiad(self):
		leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 5)
		rightA = MockSplitRead("key1", "L", "rightA", "rightFormattedRead", 42)		
		groups = {'1':([leftA],[rightA])}
		
		actual = build_pairs_from_groups(groups, MockLogger())	
		
		expected=[(leftA, rightA, 5)]
		self.assertEqual(1, len(actual))
		self.assertEqual(expected, actual)


	def test_build_pairs_from_groups(self):
		leftA = MockSplitRead("key1", "L", "leftA", "leftFormattedRead", 5)
		leftB = MockSplitRead("key1", "L", "leftB", "leftFormattedRead", 10)
		rightA = MockSplitRead("key1", "R", "rightA", "rightFormattedRead", 42)		
		rightB = MockSplitRead("key1", "R", "rightB", "rightFormattedRead", 42)		

		groups = {'1':([leftA, leftB], [rightA, rightB])}
		
		actual = build_pairs_from_groups(groups, MockLogger())	
	
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
	def __init__(self, key, side, name = "name", format = "format", distance=42):
		self._key = key
		self.name = name
		self.side = side
		self._format = format
		self._distance = distance

	def key(self):
		return self._key

	def format(self, delimiter="\t"):
		return self._format
		
	def distance(self, distance):
		return self._distance

	def __repr__(self):
		return self._name


class MockSplitReadBuilder():
	def __init__(self, d):
		self._d = d
	
	def build(self, line):
		return self._d[line]


def initParams(updates):
	params = {'read_len':30, 'name':"name", 'side':"L", 'split_len':10, 'strand':"strand", 'chr':"chr", 'position':100, 'matches':5} 
	params.update(updates);
	return params



if __name__ == "__main__":
    unittest.main() 
