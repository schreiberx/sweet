


class InfoError:

	def __init__(self, i_prefix):
		InfoError.setup(self, i_prefix)

	def setup(self, i_prefix):
		self.prefix = i_prefix

	def info(self, i_str):
		o = "INFO ["+self.prefix+"] "+i_str
		print(o)

	def error(self, i_str):
		o = "ERROR ["+self.prefix+"] "+i_str
		print(o)
		raise Exception(o)

	def hline(self):
		print("*"*80)
