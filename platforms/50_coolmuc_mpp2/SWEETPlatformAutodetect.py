
import platform

def autodetect():
	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	if platform.node()[:10] == 'mpp2-login':
		return True

	return False

if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))

