#! /usr/bin/env python3


def exec_program(progparams):

	from subprocess import Popen, PIPE
	p = Popen(progparams, stdout=PIPE, stderr=PIPE)
	output, error = p.communicate()

	error = error.decode()
	output = output.decode()

	if error != '':
		output += "\n"+error

	return output, p.returncode



if __name__ == "__main__":
	import sys
	output, code = exec_program(" ".join(sys.argv[1:]))
	print("Return code: "+str(code))
	print("Output:")
	print(output)
