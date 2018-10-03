
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

	# Fallback cheyenne platform
	return True



if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))


