
def autodetect():

    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    # Always true, since this is the fallback solution
    # This helps new user for a quick start
    return True

    # Return false to avoid any accidental errors
    #return True
    return False

if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

