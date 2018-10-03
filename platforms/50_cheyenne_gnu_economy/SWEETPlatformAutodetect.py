
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
	sweet_src_dirname = dirs[len(dirs)-4]
	print(sweet_src_dirname)

	# Autodetect based on source folder name for SWEET source
	# This helps to utilize different versions of SWEET on cheyenne
	if sweet_src_dirname=="sweet_gnu_economy":
		return True



if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))


