import platform
import socket
import sys
import os

from mule.JobGeneration import *
from mule.JobPlatformResources import *
from . import JobPlatformAutodetect

# Underscore defines symbols to be private
_job_id = None



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

    return "linuxcluster_intel"



def get_platform_resources():
    """
    Return information about hardware
    """

    h = JobPlatformResources()

    h.num_cores_per_node = 28
    # Number of nodes per job are limited
    h.num_nodes = 384
    #h.num_nodes = 60
    h.num_cores_per_socket = 14
    h.max_wallclock_seconds = 48*60*60
    return h



def jobscript_setup(jg : JobGeneration):
    """
    Setup data to generate job script
    """

    global _job_id
    _job_id = jg.runtime.getUniqueID(jg.compile, jg.unique_id_filter)
    return



def jobscript_get_header(jg : JobGeneration):
    """
    These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

    Returns
    -------
    string
    	multiline text for scripts
    """
    global _job_id

    p = jg.parallelization
    time_str = p.get_max_wallclock_seconds_hh_mm_ss()

    #
    # See https://www.lrz.de/services/compute/linux-cluster/batch_parallel/example_jobs/
    #

    content = """#! /bin/bash
#SBATCH -o """+jg.p_job_stdout_filepath+"""
#SBATCH -e """+jg.p_job_stderr_filepath+"""
#SBATCH -D """+jg.p_job_dirpath+"""
#SBATCH -J """+_job_id+"""
"""

    if p.num_nodes <= 2:
        content += "#SBATCH --clusters=cm2_tiny\n"
    else:
        content += "#SBATCH --clusters=cm2\n"

    content += """#SBATCH --get-user-env 
#SBATCH --nodes="""+str(p.num_nodes)+"""
#SBATCH --ntasks-per-node="""+str(p.num_ranks_per_node)+"""
# the above is a good match for the
# CooLMUC2 architecture.
#SBATCH --export=NONE 
#SBATCH --time="""+time_str+"""
"""

    if p.num_nodes <= 2:
    	content += "#SBATCH --partition=cm2_tiny\n"

    elif p.num_nodes <= 24:
    	content += "#SBATCH --partition=cm2_std\n"
    	content += "#SBATCH --qos=cm2_std\n"

    else:
    	content += "#SBATCH --partition=cm2_large\n"
    	content += "#SBATCH --qos=cm2_large\n"

    mule_user_email = os.getenv('MULE_USER_EMAIL')
    if mule_user_email:
        content += "#SBATCH --mail-type=BEGIN,END\n" 
        content += "#SBATCH --mail-user="+mule_user_email
    else:
        raise Exception("User email for the slurm script has not been specified. Please include it in env_vars.sh!")


    content += "\n"
    content += "module load slurm_setup\n"
 

    if False:
    	if p.force_turbo_off:
    		content += """# Try to avoid slowing down some CPUs
#SBATCH --ear=off
"""

    content += """
source /etc/profile.d/modules.sh

"""

    if jg.compile.threading != 'off':
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

    mpiexec = ''

    #
    # Only use MPI exec if we are allowed to do so
    # We shouldn't use mpiexec for validation scripts
    #
    if not p.mpiexec_disabled:
    	mpiexec = "mpiexec -n "+str(p.num_ranks)+" --perhost "+str(p.num_ranks_per_node)


    sweet_ld_library_path = os.getenv('MULE_LD_LIBRARY_PATH')
    if sweet_ld_library_path == None:
        raise Exception("Environment variable MULE_LD_LIBRARY_PATH not found!")


    content = """

# Output MPI version
echo "**************************************************"
echo "MPI Information"
echo "**************************************************"
echo "mpiexec --version"
mpiexec --version 2>&1
echo "**************************************************"


# List loaded modules
echo "**************************************************"
echo "Loaded modules"
echo "**************************************************"
echo "module list"
module list 2>&1
echo "**************************************************"

# Make sure that MULE library path is really known
export LD_LIBRARY_PATH=\""""+sweet_ld_library_path+""":$LD_LIBRARY_PATH\"
# mpiexec ... would be here without a line break
EXEC=\""""+jg.compile.getProgramPath()+"""\"
PARAMS=\""""+jg.runtime.getRuntimeOptions()+"""\"
echo \"${EXEC} ${PARAMS}\"

"""+mpiexec+""" $EXEC $PARAMS || exit 1

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

    content = f"""

SCONS="scons {jg.compile.getSConsParams()}"
echo "$SCONS"
$SCONS || exit 1
"""

    return content

