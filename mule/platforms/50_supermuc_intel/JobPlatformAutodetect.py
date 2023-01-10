
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
    if not ".sng.lrz.de" in fqdn:
    	return False

    return True


if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

