

import platform
import socket
import sys


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

	job_id = jobgeneration.runtime.getUniqueID(jobgeneration.compile)

	content = """#! /usr/bin/bash
#
# Martinium Header (job_script_header)
#
# Job id: """+job_id+"""
#

"""

	return content



#
# MPI execution command (typically mpirun)
#
def job_script_mpiexec(jobgeneration):
	return ""


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


