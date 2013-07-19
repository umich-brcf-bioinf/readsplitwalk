import unittest
import tempfile
import os

from bin.split_read import FQStanza, build_splits, write_stanzas, stanza_generator


class FQStanzaTest(unittest.TestCase):
	
	def test_parse(self):
		stanza_string = "HWI-D00196:21:H0DC6ADXX:1:1101:1306:2089 1:N:0:TAGCTT\nGCACGGTTCTGTAGTCTNCAGAAGTATCNGATNNG\n+HWI-EAS159:6:1:6:610#0/1\naaa^W\]a]^Z^]P[_[DZWR[\\\BBBBBBBBBB\n"
		actual_stanza = FQStanza.parse(stanza_string)
		
		self.assertEqual("HWI-D00196:21:H0DC6ADXX:1:1101:1306:2089_1:N:0:TAGCTT", actual_stanza.main_header)
		self.assertEqual("GCACGGTTCTGTAGTCTNCAGAAGTATCNGATNNG", actual_stanza.seq)
		self.assertEqual("+HWI-EAS159:6:1:6:610#0/1", actual_stanza.score_header)
		self.assertEqual("aaa^W\]a]^Z^]P[_[DZWR[\\\BBBBBBBBBB", actual_stanza.score)
	
	def test_split(self):
		inputHeader = "header"
		seq = "ABCDE"
		score_header = "score_header"
		score="12345"
		
		inputStanza = FQStanza(inputHeader, seq, score_header, score)
		split_position = 2
		
		(actualLeftStanza, actualRightStanza) = inputStanza.split(split_position)

		self.assertEquals("header-L-2", actualLeftStanza.main_header)
		self.assertEquals("AB", actualLeftStanza.seq)
		self.assertEquals("score_header-L-2", actualLeftStanza.score_header)
		self.assertEquals("12", actualLeftStanza.score)
		
		self.assertEquals("header-R-3", actualRightStanza.main_header)
		self.assertEquals("CDE", actualRightStanza.seq)
		self.assertEquals("score_header-R-3", actualRightStanza.score_header)
		self.assertEquals("345", actualRightStanza.score)
		
	def test_repr(self):
		input_header = "header"
		seq = "ABCDE"
		score_header = "score_header"
		score="12345"
		
		input_stanza = FQStanza(input_header, seq, score_header, score)
		
		self.assertEquals("{0}\n{1}\n{2}\n{3}".format(input_header, seq, score_header, score), str(input_stanza))
		
class SplitReadTest(unittest.TestCase):
		
	def test_build_splits(self):
		inputHeader = "header"
		seq = "ABCDE"
		score_header = "score_header"
		score="12345"
		
		inputStanza = FQStanza(inputHeader, seq, score_header, score)
		split_margin = 2
		
		stanza_list = build_splits(inputStanza, split_margin)
		
		self.assertEquals(4, len(stanza_list))
		self.assertEquals("AB", stanza_list[0].seq)
		self.assertEquals("CDE", stanza_list[1].seq)
		self.assertEquals("ABC", stanza_list[2].seq)
		self.assertEquals("DE", stanza_list[3].seq)
		
	def test_write_stanzas(self):
		stanzaA = init_stanza("headerA")
		stanzaB = init_stanza("headerB")
		writer = MockWriter()
		split_margin = 2

		write_stanzas([stanzaA, stanzaB], writer, split_margin)		

		lines = writer.lines()
		self.assertEquals(8 * 4, len(lines))
		stanzaA_headers = [s for s in lines if s.startswith("headerA")]
		self.assertEquals(4, len(stanzaA_headers))

	def test_stanza_generator(self):
		reader = MockReader(\
"""file_header1
file_header2
@h1
ABC
h1
123
@h2
DEF
+h2
456""")

		stanzas = [stanza for stanza in stanza_generator(reader, "@")]

		self.assertEquals(2, len(stanzas))
		self.assertEquals("@h1", stanzas[0].main_header)
		self.assertEquals("@h2", stanzas[1].main_header)



	def test_stanza_generator_tricky_stanza(self):
		reader = MockReader(\
"""file_header1
file_header2
@h1
ABC
@d1
123
@h2
DEF
+h2
456""")

		stanzas = [stanza for stanza in stanza_generator(reader, "@")]

		self.assertEquals(2, len(stanzas))
		self.assertEquals("@h1", stanzas[0].main_header)
		self.assertEquals("@h2", stanzas[1].main_header)


def init_stanza(header):
	return FQStanza(header, "ABCDE", "score_header", "score")
		
		
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

class MockReader():
		def __init__(self, content):
				lines = [line + "\n" for line in content.split("\n") if line != ""]
				self._iter = lines.__iter__()
				self.wasClosed = False
		
		def __iter__(self):
				return self._iter
		
		def close(self):
				self.wasClosed=True


if __name__ == "__main__":
	unittest.main()
