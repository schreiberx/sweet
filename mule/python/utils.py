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


def exec_program(progparams, shell=False, catch_output=True):

    if catch_output:
        import subprocess
        from subprocess import Popen, PIPE
        if False:
            p = Popen(progparams, stdout=PIPE, stderr=PIPE, shell=shell)
        else:
            p = subprocess.run(progparams, stdout=PIPE, stderr=PIPE, shell=shell)
        output, error = p.communicate()

        error = error.decode()
        output = output.decode()

        if error != '':
            output += "\n"+error

        return output, p.returncode

    else:
        import subprocess
        retval = subprocess.run(progparams)
        return retval.returncode


def exec_command(command):
    process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    # combine stdout and stderr
    out = out+err
    out = out.decode("utf-8")
    out = out.replace("\r", "")
    return out


def remove_file_ending(i_str):
    import os
    return os.path.splitext(i_str)[0]


