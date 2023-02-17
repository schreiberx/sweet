
def autodetect():

    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    import platform
    if platform.node() == "time-x":
        return True

    return False

if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

