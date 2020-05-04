
import platform
import os

def autodetect():

    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    if platform.node()[:8] == 'cm2login':
    	return True

    return False


if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

