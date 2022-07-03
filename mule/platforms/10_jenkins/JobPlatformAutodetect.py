import platform
import os


def autodetect():

    """
    Returns
    -------
    bool
        True if current platform matches, otherwise False
    """

    if 'CI' in os.environ:
        if os.environ['CI'] != 'true':
            return False

    if not 'JENKINS_HOME' in os.environ:
        return False

    return True

if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

