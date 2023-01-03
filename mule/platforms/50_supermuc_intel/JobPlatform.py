"""
Platform job scripts for SuperMUC, see also

https://doku.lrz.de/display/PUBLIC/Job+Processing+with+SLURM+on+SuperMUC-NG
"""

import platform
import socket
import sys
import os

from mule_local.JobGeneration import *
from mule.JobPlatformResources import *
from . import JobPlatformAutodetect

# Underscore defines symbols to be private
#_job_id = None



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

    return "supermuc_intel"



def get_platform_resources():
    """
    Return information about hardware
    """

    h = JobPlatformResources()

    """
    https://doku.lrz.de/display/PUBLIC/Hardware+of+SuperMUC-NG
    """

    h.num_cores_per_node = 48
    # Number of nodes per job are limited
    h.num_nodes = 6336
    #h.num_nodes = 60
    h.num_cores_per_socket = 24
    h.max_wallclock_seconds = 2*24*60*60
    return h



def jobscript_setup(jg : JobGeneration):
    """
    Setup data to generate job script
    """

    pass
    #_job_id = jg.runtime.getUniqueID(jg.compile, jg.unique_id_filter)



def jobscript_get_header(jg : JobGeneration):
    """
    These headers typically contain the information on e.g. Job exection, number of compute nodes, etc.

    Returns
    -------
    string
    	multiline text for scripts
    """
    #global _job_id
    _job_id = jg.runtime.getUniqueID(jg.compile, jg.unique_id_filter)

    p = jg.parallelization
    time_str = p.get_max_wallclock_seconds_hh_mm_ss()

    #
    # Next, we load the supermuc-specific configuration file which we assume to be in the HOME folder of the user
    #
    config_file = os.environ["HOME"]+"/sweet_supermuc_per_user_config.py"
    config_per_user = {}

    # Special override for unit test
    if "MULE_TEST_PLATFORMS" in os.environ:
        config_per_user['project_id'] = 'dummy_project_id'
        config_per_user['user_email'] = 'dummy_email'
    else:
        try:
            with open(config_file, "rb") as f:
                code = compile(f.read(), config_file, "exec")
                exec(code, config_per_user)

        except:
            print("*"*80)
            print("ERROR: Failed to parse '"+supermuc_config+"'")
            print("*"*80)
            print("""
Make sure that it's in the following format:

user_email = 'someemail@address.com'
project_id = 'pro123abc'

""")
            raise Exception("ERROR - stopping here")

    print("Project ID: "+config_per_user['project_id'])
    print("User Email: "+config_per_user['user_email'])

    #
    # https://doku.lrz.de/display/PUBLIC/Hardware+of+SuperMUC-NG
    #
    content = f"""#! /bin/bash
#SBATCH -J {_job_id}
#SBATCH -o {jg.p_job_stdout_filepath}
#SBATCH -e {jg.p_job_stderr_filepath}
#SBATCH -D {jg.p_job_dirpath}
"""

    """
    https://doku.lrz.de/display/PUBLIC/Access+and+Overview+of+HPC+Systems
    micro, general, large
    """

    if p.num_nodes <= 16:
        content += "#SBATCH --partition=micro\n"
        if p.max_wallclock_seconds > 48*60*60:
            raise Exception("Max. wallclock time exceeds maximum time")

    elif p.num_nodes <= 768:
        content += "#SBATCH --partition=general\n"

        if p.max_wallclock_seconds > 48*60*60:
            raise Exception("Max. wallclock time exceeds maximum time")
    else:
        content += "#SBATCH --partition=large\n"

        if p.max_wallclock_seconds > 24*60*60:
            raise Exception("Max. wallclock time exceeds maximum time")

    content += f"""
#SBATCH --nodes={p.num_nodes}
#SBATCH --ntasks-per-node={p.num_ranks_per_node}
#SBATCH --mail-type=end 
#SBATCH --mail-user={config_per_user['user_email']}
#SBATCH --account={config_per_user['project_id']}
#SBATCH --time={time_str}
#SBATCH --no-requeue
#SBATCH --get-user-env 
#SBATCH --export=NONE 
"""

    if True:
    	if p.force_turbo_off:
    		content += """# Try to avoid slowing down some CPUs
#SBATCH --ear=off
"""

    content += "\n"
    content += "module load slurm_setup\n"
 


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

    content = """

SCONS="scons """+jg.compile.getSConsParams()+' -j 4"'+"""
echo "$SCONS"
$SCONS || exit 1
"""

    return content

