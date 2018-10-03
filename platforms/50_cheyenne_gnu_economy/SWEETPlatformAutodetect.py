
import platform
import socket

def autodetect():
	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	fqdn = socket.getfqdn()
	return ".cheyenne" in fqdn


if __name__ == "__main__":
	print("Autodetect: "+str(autodetect()))


