
import platform
import os

def autodetect():

    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    n = platform.node()

    # Login node
    if n == "vmschulz8":
        return True

    if len(n) == 5:
        if n[:3] in ['odr', 'rpi']:
            return True

    return False


if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

