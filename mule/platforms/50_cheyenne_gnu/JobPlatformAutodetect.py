
import platform
import socket
import os

def autodetect():
	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	fqdn = socket.getfqdn()
	if not ".cheyenne" in fqdn:
		return False

	dirs = os.path.abspath(__file__).split('/')
	sweet_src_dirname = dirs[len(dirs)-5]

	# Autodetect based on source folder name for MULE source
	# This helps to utilize different versions of MULE on cheyenne
	if sweet_src_dirname=="sweet_gnu":
		return True

	return False


if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))


