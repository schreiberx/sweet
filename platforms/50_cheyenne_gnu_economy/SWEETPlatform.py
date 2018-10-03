import platform
import socket
import sys

from SWEETPlatformResources import *
from SWEETJobGeneration import *
from . import SWEETPlatformAutodetect

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



def p_gen_script_info(jobgeneration : SWEETJobGeneration):
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

	return SWEETPlatformAutodetect.get_platform_autodetect()



def get_platform_id():
	"""
	Return platform ID

	Returns
	-------
	string
		unique ID of platform
	"""

	return "cheyenne_gnu_economy"



def get_platform_resources():
	"""
	Return information about hardware
	"""

	h = SWEETPlatformResources()

	h.num_cores_per_node = 36
	# Physical number of nodes
	h.num_nodes = 4032
	h.num_cores_per_socket = 18
	# 12h limit
	h.max_wallclock_seconds = 60*60*12
	return h



def jobscript_setup(jobgeneration : SWEETJobGeneration):
	"""
	Setup data to generate job script
	"""

	global _job_id
	_job_id = jobgeneration.runtime.getUniqueID(jobgeneration.compile)
	return



def jobscript_get_header(jobgeneration : SWEETJobGeneration):
	"""
	These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

	Returns
	-------
	string
		multiline text for scripts
	"""
	global _job_id

	p = jobgeneration.parallelization
	h = jobgeneration.platform_hardware

	def get_time_str(secs):
		# seconds
		s = int(secs)
		m = s // 60
		s = s % 60
		h = m // 60
		m = m % 60

		stest = h*60*60 + m*60 + s
		if int(secs) != stest:
			print(secs)
			print(stest)
			raise Exception("Internal error!")

		time_str = str(h).zfill(2)+":"+str(m).zfill(2)+":"+str(s).zfill(2)

	time_str = get_time_str(p.max_wallclock_seconds)
	
	#
	# See https://www.lrz.de/services/compute/linux-cluster/batch_parallel/example_jobs/
	#
	content = """#! /bin/bash
#
## project code
#PBS -A NCIS0002
## economy queue
#PBS -q economy
## wall-clock time (hrs:mins:secs)
#PBS -l walltime="""+time_str+"""
## select: number of nodes
## ncpus: number of CPUs per node
## mpiprocs: number of ranks per node
#PBS -l select="""+str(p.num_nodes)+""":ncpus="""+str(h.num_cores_per_node)+""":mpiprocs="""+str(p.num_ranks_per_node)+""":ompthreads="""+str(p.threads_per_rank)+"\n"

	#"default": 2301000 
	#"turbo": 2301000
	#"rated": 2300000
	#"slow": 1200000
	if jobgeneration.cluster.force_turbo_off:
		content += "#PBS -l select=cpufreq=2300000\n"

	content += """#
#PBS -N """+_job_id[0:100]+"""
#PBS -o """+jobgeneration.p_job_stdout_filepath+"""
#PBS -o """+jobgeneration.p_job_stderr_filepath+"""

source /etc/profile.d/modules.sh

export OMP_NUM_THREADS="""+str(p.num_threads_per_rank)+"""

module load openmpi
"""+("module load mkl" if mkl else "")+"""


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

	mpiexec = ''

	#
	# Only use MPI exec if we are allowed to do so
	# We shouldn't use mpiexec for validation scripts
	#
	if not p.mpiexec_disabled:
		mpiexec = "mpiexec_mpt -n "+str(p.num_ranks)

		mpiexec += " omplace "
		mpiexec += " -nt "+str(p.threads_per_rank)+" "
		mpiexec += " -vv "


	content = """

"""+p_gen_script_info(jobgeneration)+"""

# mpiexec ... would be here without a line break
EXEC=\"$SWEET_ROOT/build/"""+jobgeneration.compile.getProgramName()+"""\"
PARAMS=\""""+jobgeneration.runtime.getRuntimeOptions()+"""\"
echo \"${EXEC} ${PARAMS}\"

"""+mpiexec+""" $EXEC $PARAMS

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
$SCONS || exit 1
"""

	return content

