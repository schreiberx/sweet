import platform
import os


def autodetect():

    """
    Returns
    -------
    bool
        True if current platform matches, otherwise False
    """

    if 'GITLAB_CI' in os.environ:
        if os.environ['GITLAB_CI'] != 'true':
            return False

    return True

if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

