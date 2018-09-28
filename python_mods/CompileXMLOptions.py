from xml.etree.ElementTree import XMLParser
import sys
import os
import commands



class CXMLParser:
	depth = 0
	compiler_scope = False

	def start(self, tag, attrib):
		if self.depth == 0:
			pass

		elif self.depth == 1:
			if tag in ['compiler', 'c3s']:
				self.compiler_scope = True

		elif self.depth == 2:
			if self.compiler_scope:
				if tag == 'param':
					if attrib['value'] != '*':
						attrib['name'] = attrib['name'].replace('-', '_')
						self.env[attrib['name']] = attrib['value']

		self.depth += 1

	def end(self, tag):
		self.depth -= 1

		if self.depth == 1:
			if tag == 'compiler':
				self.compiler_scope = False

	def data(self, data):
		pass

	def close(self):
		pass

	def __init__(self, i_env):
		self.env=i_env



def load(file, env):
	cXMLParser = CXMLParser(env)

	parser = XMLParser(target=cXMLParser)

	f = open(file, 'r+')
	read_data = f.read();
	f.close()

	parser.feed(read_data)
	parser.close()
	
	return cXMLParser.env
