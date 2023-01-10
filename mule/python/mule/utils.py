#! /usr/bin/env python3

def hostname():

    # First, try with environment variables
    import os
    if 'HOSTNAME' in os.environ:
        hostname = os.environ['HOSTNAME']

        if hostname != "":
            return hostname

    # Second, try with platform node
    import platform
    hostname = platform.node()

    if hostname is None:
        raise Exception("Can't detect hostname")

    return hostname
