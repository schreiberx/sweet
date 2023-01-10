import platform
import socket
import sys
import os

from mule_local.JobGeneration import *
from mule.JobPlatformResources import *
from . import JobPlatformAutodetect

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
    return """#
# Generating function: """+_whoami(2)+"""
# Platform: """+get_platform_id()+"""
# Job id: """+jg.getUniqueID()+"""
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

    return "cheyenne_gcc"



def get_platform_resources():
    """
    Return information about hardware
    """

    r = JobPlatformResources()

    r.num_cores_per_node = 36

    # Physical number of nodes, maybe the limit is different
    r.num_nodes = 4032

    r.num_cores_per_socket = 18

    # 12h limit
    r.max_wallclock_seconds = 60*60*12
    return r



def jobscript_setup(jg : JobGeneration):
    """
    Setup data to generate job script
    """

    return



def jobscript_get_header(jg : JobGeneration):
    """
    These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

    Returns
    -------
    string
    	multiline text for scripts
    """
    job_id = jg.getUniqueID()

    p = jg.parallelization

    time_str = p.get_max_wallclock_seconds_hh_mm_ss()
    
    # Available queues:
    # premium	(only use this in extreme cases)
    # regular
    # economy
    queue = 'economy'

    # Use regular queue if we need more than 32 nodes
    # Otherwise, the job doesn't seem to be scheduled

    if p.num_nodes >= 32:
    	queue = 'premium'
    elif p.num_nodes >= 16:
    	queue = 'regular'

    #
    # See https://www.lrz.de/services/compute/linux-cluster/batch_parallel/example_jobs/
    #
    content = """#! /bin/bash
#
## project code
#PBS -A NCIS0002
#PBS -q """+queue+"""
## wall-clock time (hrs:mins:secs)
#PBS -l walltime="""+time_str+"""
## select: number of nodes
## ncpus: number of CPUs per node
## mpiprocs: number of ranks per node
#PBS -l select="""+str(p.num_nodes)+""":ncpus="""+str(p.num_cores_per_node)+""":mpiprocs="""+str(p.num_ranks_per_node)+""":ompthreads="""+str(p.num_threads_per_rank)+"\n"

    #"default": 2301000 
    #"turbo": 2301000
    #"rated": 2300000
    #"slow": 1200000
    if p.force_turbo_off:
    	content += "#PBS -l select=cpufreq=2300000\n"

    content += """#
#PBS -N """+job_id[0:100]+"""
#PBS -o """+jg.p_job_stdout_filepath+"""
#PBS -e """+jg.p_job_stderr_filepath+"""

#source /etc/profile.d/modules.sh

#module load openmpi
"""+("module load mkl" if jg.compile.mkl==True or jg.compile.mkl=='enable' else "")+"""


"""+p_gen_script_info(jg)+"""


echo
echo "hostname"
hostname
echo

echo
echo "lscpu -e"
lscpu -e 
echo

echo
echo "CPU Frequencies (uniquely reduced):"
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq | sort -u
echo

"""

    if jg.compile.threading != 'off':
    	content += """
export OMP_NUM_THREADS="""+str(p.num_threads_per_rank)+"""
export OMP_DISPLAY_ENV=VERBOSE
"""

    if p.core_oversubscription:
    	raise Exception("Not supported with this script!")
    else:
    	if p.core_affinity != None:
    		content += "\necho \"Affnity: "+str(p.core_affinity)+"\"\n"
    		if p.core_affinity == 'compact':
    			content += "source $MULE_ROOT/platforms/bin/setup_omp_places.sh nooversubscription close\n"
    			#content += "\nexport OMP_PROC_BIND=close\n"
    		elif p.core_affinity == 'scatter':
    			raise Exception("Affinity '"+str(p.core_affinity)+"' not supported")
    			content += "\nexport OMP_PROC_BIND=spread\n"
    		else:
    			raise Exception("Affinity '"+str(p.core_affinity)+"' not supported")

    		content += "\n"

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

    content += """

EXEC=\""""+jg.compile.getProgramPath()+"""\"
PARAMS=\""""+jg.runtime.getRuntimeOptions()+"""\"

"""

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

    mpiexec = ""

    #
    # Only use MPI exec if we are allowed to do so
    # We shouldn't use mpiexec for validation scripts
    #
    if not p.mpiexec_disabled:
    	# Use mpiexec_mpt for Intel MPI
    	#mpiexec = "mpiexec_mpt -n "+str(p.num_ranks)

    	# Use mpiexec for GNU
    	if jg.compile.sweet_mpi == 'enable':
    		mpiexec = "mpiexec_mpt -n "+str(p.num_ranks)

    		mpiexec += " omplace "
    		#mpiexec += " -nt "+str(p.num_threads_per_rank)+" "
    		mpiexec += " -nt "+str(p.num_cores_per_rank)+" "
    		# Don't know if intel mode really works with gcc
    		mpiexec += " -tm intel "
    		mpiexec += " -vv"
    		if mpiexec[-1] != ' ':
    			mpiexec += ' '

    #
    # Fix the mess on Cheyenne!
    #
    # We prefix the current LD_LIBRARY_PATH with the one from the shell where the job was submitted
    # This is required since Cheyenne scripts mess around with the existing path in a way
    # which results in e.g. the system-wide installed fftw to be loaded.
    #
    # What we basically accomplish here is to suggest to really first
    # lookup the MULE local_software/local/lib directory, then the system libraries
    #
    sweet_ld_library_path = os.getenv('MULE_LD_LIBRARY_PATH')

    if sweet_ld_library_path == None:
    	raise Exception("Environment variable MULE_LD_LIBRARY_PATH not found!")

    content = """


# Make sure that MULE library path is really known
export LD_LIBRARY_PATH=\""""+sweet_ld_library_path+""":$LD_LIBRARY_PATH\"


echo
echo "LD_LIBRARY_PATH"
echo "${LD_LIBRARY_PATH}"
echo

echo
echo "ldd"
ldd $EXEC
echo

E=\""""+mpiexec+"""${EXEC} ${PARAMS}\"

echo
echo "Executing..."
echo "$E"
$E || exit 1

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

    content = """
echo
echo "CPU Frequencies (uniquely reduced):"
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq | sort -u
echo

"""

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
    content = ""

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

