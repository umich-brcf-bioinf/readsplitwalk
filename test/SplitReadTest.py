#! /usr/bin/env python

import sys ; sys.path.insert(0, "../bin")
#input example: /ccmb/BioinfCore_ARCHIVE/Vol01/ccmbCoreProjects_Archieve_YBAI/Run_46/Probe/s6.unmapped.fa
#output example: /ccmb/BioinfCore_ARCHIVE/Vol01/ccmbCoreProjects_Archieve_YBAI/Run_46/Probe/s6.unmapped.fa.split.left15_right17

import unittest
import tempfile
import os

from SplitRead import FQStanza


class FQStanzaTest(unittest.TestCase):
	
	def test_init(self):
		stanza_string = "@HWI-EAS159:6:1:6:610#0/1\nGCACGGTTCTGTAGTCTNCAGAAGTATCNGATNNG\n+HWI-EAS159:6:1:6:610#0/1\naaa^W\]a]^Z^]P[_[DZWR[\\\BBBBBBBBBB\n"
		actual_stanza = FQStanza(stanza_string)
		
		self.assertEqual("@HWI-EAS159:6:1:6:610#0/1", actual_stanza.main_header)
		self.assertEqual("GCACGGTTCTGTAGTCTNCAGAAGTATCNGATNNG", actual_stanza.seq)
		self.assertEqual("+HWI-EAS159:6:1:6:610#0/1", actual_stanza.score_header)
		self.assertEqual("aaa^W\]a]^Z^]P[_[DZWR[\\\BBBBBBBBBB", actual_stanza.score)



if __name__ == "__main__":
	unittest.main() 
	print "done."
