
import platform

def autodetect():

	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	prefix = "pedrosp-envy"
	if platform.node()[:len(prefix)] == prefix:
		return True

	return False


if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))

