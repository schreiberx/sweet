

import platform
import socket



#
#
#
def platform_autodetect():
	if platform.node() == 'martinium':
		return True
	else:
		return False


#
# Return an identifier string if system detected, otherwise 'None'
#
def platform_id():
	return 'martinium'



#
# Header for scripts
# These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.
#
def job_script_header(jobgeneration, dirname):
	content = ""

	content = """#! /usr/bin/bash
#
# Yellowstone Header, using external compile scripts
#

"""
	f = open('compile_yellowstone.sh', 'w')
	f.write("#! /bin/bash\n")
	f.write("\n")

	if jobgeneration.compile.mkl == 'enable':
		f.write("module load mkl\n")
		f.write("\n")

	f.write("scons "+jobgeneration.compile.getSConsParams()+'\n')
	f.write("\n")
	os.chmod(fn, 0o755)

	return content


#
# Prefix before the executable is executed
#
def job_script_exec_prefix(jobgeneration):
	content = ""

	content += """
#
# Martinium compilation (job_script_exec_prefix)
#
SCONS="scons """+jobgeneration.compile.getSConsParams()+' -j 4"'+"""
echo "$SCONS"
$SCONS
"""

	return content


