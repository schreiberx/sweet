import platform
import socket
import sys

from SWEETPlatformResources import *
from SWEETJobGeneration import *
from . import SWEETPlatformAutodetect

# Underscore defines symbols to be private
_job_id = None



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

	return "coolmuc_mpp2_gnu"



def get_platform_resources():
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



def jobscript_setup(j : SWEETJobGeneration):
	"""
	Setup data to generate job script
	"""

	global _job_id
	_job_id = j.runtime.getUniqueID(j.compile)
	return



def jobscript_get_header(j : SWEETJobGeneration):
	"""
	These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

	Returns
	-------
	string
		multiline text for scripts
	"""
	global _job_id

	p = j.parallelization
	time_str = p.get_max_wallclock_seconds_hh_mm_ss()

	#
	# See https://www.lrz.de/services/compute/linux-cluster/batch_parallel/example_jobs/
	#
	content = """#! /bin/bash
#SBATCH -o """+j.p_job_stdout_filepath+"""
#SBATCH -e """+j.p_job_stderr_filepath+"""
#SBATCH -D """+j.p_job_dirpath+"""
#SBATCH -J """+_job_id+"""
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
"""

	if False:
		if p.force_turbo_off:
			content += """# Try to avoid slowing down CPUs
#SBATCH --cpu-freq=Performance
"""

	content += """
source /etc/profile.d/modules.sh

"""

	if j.compile.threading != 'off':
		content += """
export OMP_NUM_THREADS="""+str(p.num_threads_per_rank)+"""
"""

	if p.core_oversubscription:
		raise Exception("Not supported with this script!")

	if p.core_affinity != None:
		
		content += "\necho \"Affnity: "+str(p.core_affinity)+"\"\n"
		if p.core_affinity == 'compact':
			content += "\nexport OMP_PROC_BIND=close\n"
		elif p.core_affinity == 'scatter':
			content += "\nexport OMP_PROC_BIND=spread\n"
		else:
			raise Exception("Affinity '"+str(p.core_affinity)+"' not supported")

	return content






def jobscript_get_exec_prefix(j : SWEETJobGeneration):
	"""
	Prefix before executable

	Returns
	-------
	string
		multiline text for scripts
	"""

	p = j.parallelization
	content = ""

	return content



def jobscript_get_exec_command(j : SWEETJobGeneration):
	"""
	Prefix to executable command

	Returns
	-------
	string
		multiline text for scripts
	"""

	p = j.parallelization

	mpiexec = ''

	#
	# Only use MPI exec if we are allowed to do so
	# We shouldn't use mpiexec for validation scripts
	#
	if not p.mpiexec_disabled:
		mpiexec = "mpiexec -n "+str(p.num_ranks)+" --perhost "+str(p.num_ranks_per_node)

	content = """

# mpiexec ... would be here without a line break
EXEC=\"$SWEET_ROOT/build/"""+j.compile.getProgramName()+"""\"
PARAMS=\""""+j.runtime.getRuntimeOptions()+"""\"
echo \"${EXEC} ${PARAMS}\"

"""+mpiexec+""" $EXEC $PARAMS

"""

	return content



def jobscript_get_exec_suffix(j : SWEETJobGeneration):
	"""
	Suffix before executable

	Returns
	-------
	string
		multiline text for scripts
	"""

	content = ""
	return content



def jobscript_get_footer(j : SWEETJobGeneration):
	"""
	Footer at very end of job script

	Returns
	-------
	string
		multiline text for scripts
	"""

	content = ""
	return content



def jobscript_get_compile_command(j : SWEETJobGeneration):
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

SCONS="scons """+j.compile.getSConsParams()+' -j 4"'+"""
echo "$SCONS"
$SCONS || exit 1
"""

	return content

