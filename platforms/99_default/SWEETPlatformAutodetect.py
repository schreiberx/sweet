
def autodetect():

	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	# Always true, since this is the fallback solution
	#return True

	# Return false to avoid any accidental errors
	return False

if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))

