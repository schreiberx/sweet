#! /usr/bin/env python3
"""
Some first string scaling tests for preliminary results
"""

import sys
from itertools import product

from mule.JobGeneration import JobGeneration
from mule.JobParallelizationDimOptions import JobParallelizationDimOptions
from mule.sdc import getSDCSetup

p = JobGeneration()
verbose = True

p.unique_id_filter = [
    'compile',
    'runtime']

# Main parameters
nPointsSpace = 256
dt = 128/nPointsSpace*300.0
nSteps = 100
nProcSpace = [1, 2, 4, 8, 16, 32]
nProcTime = [1, 4]
spaceTimePar = {
    nT: [nS for nS in nProcSpace if nS*nT <= max(nProcSpace)]
    for nT in nProcTime
}

# Compilation settings
p.compile.program = 'programs/pde_sweSphere'
p.compile.mode = 'release'
p.compile.lapack = 'enable'
p.compile.mkl = 'disable'

p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'

p.compile.threading = 'omp'
p.compile.rexi_thread_parallel_sum = 'disable'
p.compile.parallel_sdc_par_model = 'omp'

p.compilecommand_in_jobscript = False

# Runtime settings
p.runtime.timestepping_method = 'ln_imex_sdc'
p.runtime.paramsSDC = getSDCSetup(
    nNodes=4,
    nIter=3,
    nodeType='RADAU-RIGHT', 
    qDeltaImplicit='BEPAR', 
    qDeltaExplicit='PIC', 
    qDeltaInitial='BEPAR',
    diagonal=True,
    initialSweepType="QDELTA",
    useEndUpdate=False
)

p.runtime.space_res_spectral = nPointsSpace
p.runtime.space_res_physical = -1

p.runtime.benchmark_name = "galewsky"
p.runtime.timestep_size = dt
p.runtime.max_simulation_time = nSteps*dt

p.runtime.verbosity = 0
p.runtime.instability_checks = 0
p.runtime.viscosity = 0.0

# Parallelization settings
p.parallelization.core_oversubscription = False
p.parallelization.core_affinity = 'compact'

if __name__ == "__main__":

    for nT, listNProcSpace in spaceTimePar.items():
        for nS in listNProcSpace:
            
            # Update TIME parallelization
            ptime = JobParallelizationDimOptions('time')
            ptime.num_cores_per_rank = nT
            ptime.num_threads_per_rank = nT
            ptime.num_ranks = 1

            pspace = JobParallelizationDimOptions('space')
            pspace.num_cores_per_rank = nS
            pspace.num_threads_per_rank = nS
            pspace.num_ranks = 1

            # Setup parallelization
            p.setup_parallelization([pspace, ptime])

            if verbose:
                pspace.print()
                p.parallelization.print()

            p.parallelization.max_wallclock_seconds = 60*60*2

            p.gen_jobscript_directory()

    p.write_compilecommands()
           

            

            

            

    
