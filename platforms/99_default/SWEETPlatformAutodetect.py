
def autodetect():

	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	# Always true, since this is the fallback solution
	return True

if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))

