#! /usr/bin/env python

import unittest
import tempfile
import os

from SplitFile import splitfile, FileSystem

#from IdentifyPairsTest import MockWriter

class SplitFileTest(unittest.TestCase):
	
	def test_splitfile_identity(self):
		col_index = 2
		output_path = "outfile/path/"
		delim = "\|"
		test_body = "read1|+|chr14\nread2|+|chr14\nread3|+|chr14"
		input_filename = "tempfile.txt"
		input_filepath = "/foo/bar/" + input_filename
		file_system = MockFileSystem({input_filepath:test_body})
		
		splitfile(file_system, input_filepath, output_path, col_index, delim)
		
		writers = file_system.writers
		self.assertEqual(1, len(writers))
		newfile_name = "{0}{1}.{2}.prt".format(output_path, input_filename, "chr14")
		self.assertEqual(test_body.splitlines(), writers[newfile_name].lines())


	def test_splitfile_twoFiles(self):
		col_index = 2
		output_path = ""
		delim = "\|"
		test_body_a = "read1|+|chr14\nread2|+|chr14\nread3|+|chr14\n"
		test_body_b = "read4|+|chr15"
		input_filepath = "tempfile.txt"	
		file_system = MockFileSystem({input_filepath : test_body_a + test_body_b})

		splitfile(file_system, input_filepath , output_path, col_index, delim)
	
		writers = file_system.writers
		self.assertEqual(2, len(writers))

		newfile_name1 = "{0}.{1}.prt".format(input_filepath, "chr14")
		newfile_name2 = "{0}.{1}.prt".format(input_filepath, "chr15")
		self.assertEqual(test_body_a.splitlines(), writers[newfile_name1].lines())
		self.assertEqual(test_body_b.splitlines(), writers[newfile_name2].lines())
		
		for key in writers:
			self.assertEqual(True, writers[key].wasClosed)
			
	def test_splitfile_funkyRE(self):
		col_index = 2
		output_path = ""
		delim = "\s+"
		test_body_a = "read1   + chr14\nread2  +   chr14\nread3      +\tchr14"
		input_filepath = "tempfile.txt"	
		file_system = MockFileSystem({input_filepath : test_body_a})

		splitfile(file_system, input_filepath , output_path, col_index, delim)
	
		writers = file_system.writers
		self.assertEqual(1, len(writers))

		newfile_name1 = "{0}.{1}.prt".format(input_filepath, "chr14")
		self.assertEqual(test_body_a.splitlines(), writers[newfile_name1].lines())
		
		for key in writers:
			self.assertEqual(True, writers[key].wasClosed)


class MockFileSystem():
	
	def __init__(self, reader_dict):
		self._readers = reader_dict
		self.writers = {}
	
	def open_file(self, filename, mode):
		if mode == "r":
			return MockReader(self._readers[filename])
		elif mode == "w":
			self.writers[filename]=MockWriter()
			return self.writers[filename]

class MockReader():
	
	def __init__(self, content):
		lines = [line + "\n" for line in content.split("\n") if line != ""]
		self._iter = lines.__iter__()
		self.wasClosed = False

	def __iter__(self):
		return self._iter

	def next(self):
		self._iter.next()

	def close(self):
		self.wasClosed=True


class MockWriter():
	def __init__(self):
		self._content = []
		self.wasClosed = False

	def write(self, content):
		self._content.append(content)
		
	def lines(self):
		return "".join(self._content).splitlines()

	def close(self):
		self.wasClosed = True

if __name__ == "__main__":
	unittest.main() 
	print "done."