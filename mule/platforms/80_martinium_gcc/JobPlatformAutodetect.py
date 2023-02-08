
import platform
import os

def autodetect():

    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    prefix = "martinium"
    if platform.node()[:len(prefix)] != prefix:
    	return False

    dirs = os.path.abspath(__file__).split('/')
    sweet_src_dirname = dirs[len(dirs)-5]

    # Autodetect based on source folder name for MULE source
    # This helps to utilize different versions of MULE on cheyenne
    if sweet_src_dirname=="sweet_gcc_economy":
    	return True


    if sweet_src_dirname=="sweet":
    	return True

#    return True	### TODO: Remove this for more fine-granular detecction
    return False


if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

