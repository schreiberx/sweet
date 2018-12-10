import platform
import socket
import sys

from mule_local.JobMule import *
from . import JobPlatformAutodetect

import multiprocessing


# Underscore defines symbols to be private
_job_id = None

def _whoami(depth=1):
	"""
	String of function name to recycle code

	https://www.oreilly.com/library/view/python-cookbook/0596001673/ch14s08.html

	Returns
	-------
	string
		Return function name
	"""
	return sys._getframe(depth).f_code.co_name



def p_gen_script_info(jg : JobGeneration):
	global _job_id

	return """#
# Generating function: """+_whoami(2)+"""
# Platform: """+get_platform_id()+"""
# Job id: """+_job_id+"""
#
"""



def get_platform_autodetect():
	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	return JobPlatformAutodetect.autodetect()



def get_platform_id():
	"""
	Return platform ID

	Returns
	-------
	string
		unique ID of platform
	"""

	return "martinium_gnu"


def get_platform_resources():
	"""
	Return information about hardware
	"""

	h = JobPlatformResources()

	#h.num_cores_per_node = multiprocessing.cpu_count()
	# There are only 2 cores, not 4 !!!
	h.num_cores_per_node = 2
	h.num_nodes = 1

	# TODO: So far, we only assume a single socket system as a fallback
	h.num_cores_per_socket = h.num_cores_per_node

	return h



def jobscript_setup(jg : JobGeneration):
	"""
	Setup data to generate job script
	"""

	global _job_id
	_job_id = jg.runtime.getUniqueID(jg.compile)
	return



def jobscript_get_header(jg : JobGeneration):
	"""
	These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

	Returns
	-------
	string
		multiline text for scripts
	"""
	p = jg.parallelization

	content = """#! /bin/bash

export OMP_NUM_THREADS="""+str(p.num_threads_per_rank)+"""

"""

	return content




def jobscript_get_exec_prefix(jg : JobGeneration):
	"""
	Prefix before executable

	Returns
	-------
	string
		multiline text for scripts
	"""
	content = ""

	content += jg.runtime.get_jobscript_plan_exec_prefix(jg.compile, jg.runtime)

	return content



def jobscript_get_exec_command(jg : JobGeneration):
	"""
	Prefix to executable command

	Returns
	-------
	string
		multiline text for scripts
	"""

	p = jg.parallelization

	limit_to_physical_cores = "$MULE_ROOT/platforms/bin/exec_on_physical_cores.sh "

	content = """

"""+p_gen_script_info(jg)+"""

# mpiexec ... would be here without a line break
EXEC=\""""+jg.compile.getProgramPath()+"""\"
PARAMS=\""""+jg.runtime.getRuntimeOptions()+"""\"
echo \"${EXEC} ${PARAMS}\"

# Ensure that program is only executed on physical cores!
"""+limit_to_physical_cores+""" $EXEC $PARAMS || exit 1

"""

	return content



def jobscript_get_exec_suffix(jg : JobGeneration):
	"""
	Suffix before executable

	Returns
	-------
	string
		multiline text for scripts
	"""

	content = ""
	content += jg.runtime.get_jobscript_plan_exec_suffix(jg.compile, jg.runtime)

	return content



def jobscript_get_footer(jg : JobGeneration):
	"""
	Footer at very end of job script

	Returns
	-------
	string
		multiline text for scripts
	"""
	content = """

"""+p_gen_script_info(jg)+"""

"""

	return content



def jobscript_get_compile_command(jg : JobGeneration):
	"""
	Compile command(s)

	This is separated here to put it either
	* into the job script (handy for workstations)
	or
	* into a separate compile file (handy for clusters)

	Returns
	-------
	string
		multiline text with compile command to generate executable
	"""

	content = """

SCONS="scons """+jg.compile.getSConsParams()+' -j 4"'+"""
echo "$SCONS"
$SCONS || exit 1
"""

	return content

