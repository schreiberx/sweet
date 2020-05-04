#! /usr/bin/env python3


def exec_program(progparams, shell=True, catch_output=True):

    if catch_output:
        from subprocess import Popen, PIPE
        p = Popen(progparams, stdout=PIPE, stderr=PIPE, shell=shell)
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


if __name__ == "__main__":
    import sys
    output, code = exec_program(" ".join(sys.argv[1:]))
    #output, code = exec_program(sys.argv[1:])
    print("Return code: "+str(code))
    print("Output:")
    print(output)
