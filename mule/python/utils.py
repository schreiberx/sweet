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


def exec_command(command):
    process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    # combine stdout and stderr
    out = out+err
    out = out.decode("utf-8")
    out = out.replace("\r", "")
    return out

