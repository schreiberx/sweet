import platform
import socket
import sys

from SWEETPlatformResources import *
import SWEETJobGeneration

job_id = None

def p_whoami(depth=1):
	"""
	String of function name to recycle code

	https://www.oreilly.com/library/view/python-cookbook/0596001673/ch14s08.html

	Returns
	-------
	string
		Return function name
	"""
	return sys._getframe(depth).f_code.co_name



def p_gen_script_info(jobgeneration : SWEETJobGeneration):
	global job_id

	return """#
# Platform: """+get_platform_id()+"""
# Generating function: """+p_whoami(2)+"""
# Job id: """+job_id+"""
#
"""



def get_platform_autodetect():
	"""
	Returns
	-------
	bool
		True if current platform matches, otherwise False
	"""

	if platform.node()[:10] == 'mpp2-login':
		return True

	return False



def get_platform_id():
	"""
	Return platform ID

	Returns
	-------
	string
		unique ID of platform
	"""

	return "coolmuc_mpp2"



def get_platform_hardware():
	"""
	Return information about hardware
	"""

	h = SWEETPlatformResources()

	h.num_cores_per_node = 28
	# Number of nodes per job are limited
	#h.num_nodes = 384
	h.num_nodes = 60
	h.num_cores_per_socket = 14
	h.max_wallclock_seconds = 48*60*60
	return h



def jobscript_setup(jobgeneration : SWEETJobGeneration):
	"""
	Setup data to generate job script
	"""

	global job_id
	job_id = jobgeneration.runtime.getUniqueID(jobgeneration.compile)
	return



def jobscript_get_header(jobgeneration : SWEETJobGeneration):
	"""
	These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

	Returns
	-------
	string
		multiline text for scripts
	"""
	global job_id

	p = jobgeneration.parallelization

	from time import gmtime, strftime
	time_str = strftime("%H:%M:%S", gmtime(p.max_wallclock_seconds))
	
	#
	# See https://www.lrz.de/services/compute/linux-cluster/batch_parallel/example_jobs/
	#
	content = """#! /bin/bash
#SBATCH -o """+jobgeneration.p_jobscript_stdout_filepath+"""/output.out
#SBATCH -D """+jobgeneration.p_jobscript_dirpath+"""
#SBATCH -J """+job_id+"""
#SBATCH --get-user-env 
#SBATCH --clusters=mpp2
#SBATCH --ntasks="""+str(p.num_ranks)+"""
#SBATCH --cpus-per-task="""+str(p.num_cores_per_rank)+"""
# the above is a good match for the
# CooLMUC2 architecture.
#SBATCH --mail-type=end 
#SBATCH --mail-user=schreiberx@gmail.com
#SBATCH --export=NONE 
#SBATCH --time="""+time_str+"""
source /etc/profile.d/modules.sh

"""+p_gen_script_info(jobgeneration)+"""

"""

	return content






def jobscript_get_exec_prefix(jobgeneration : SWEETJobGeneration):
	"""
	Prefix before executable

	Returns
	-------
	string
		multiline text for scripts
	"""

	p = jobgeneration.parallelization

	content = """

"""+p_gen_script_info(jobgeneration)+"""

export OMP_NUM_THREADS="""+str(p.num_threads_per_rank)+"""

"""

	return content



def jobscript_get_exec_command(jobgeneration : SWEETJobGeneration):
	"""
	Prefix to executable command

	Returns
	-------
	string
		multiline text for scripts
	"""

	p = jobgeneration.parallelization

	content = """

"""+p_gen_script_info(jobgeneration)+"""

# mpiexec ... would be here without a line break
EXEC=\"$SWEETROOT/build/"""+jobgeneration.compile.getProgramName()+""" """+jobgeneration.runtime.getRuntimeOptions()+"""\"
echo \"$EXEC\"

mpiexec -n """+str(p.num_ranks)+""" --perhost """+str(p.num_ranks_per_node)+""" "$EXEC"

"""

	return content



def jobscript_get_exec_suffix(jobgeneration : SWEETJobGeneration):
	"""
	Suffix before executable

	Returns
	-------
	string
		multiline text for scripts
	"""

	content = """

"""+p_gen_script_info(jobgeneration)+"""

"""

	return content



def jobscript_get_footer(jobgeneration : SWEETJobGeneration):
	"""
	Footer at very end of job script

	Returns
	-------
	string
		multiline text for scripts
	"""
	content = """

"""+p_gen_script_info(jobgeneration)+"""

"""

	return content



def jobscript_get_compile_command(jobgeneration : SWEETJobGeneration):
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

SCONS="scons """+jobgeneration.compile.getSConsParams()+' -j 4"'+"""
echo "$SCONS"
$SCONS
"""

	return content

