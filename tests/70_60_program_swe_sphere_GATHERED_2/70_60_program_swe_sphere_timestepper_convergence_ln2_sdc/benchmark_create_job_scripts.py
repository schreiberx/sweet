#! /usr/bin/env python3
from mule.JobMule import JobGeneration
from mule.sdc import getSDCSetup

jg = JobGeneration()


#
# Run simulation on plane or sphere
#
jg.compile.program = 'programs/pde_sweSphere'

jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'

# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = 64
#jg.runtime.space_res_spectral = 128
jg.runtime.space_res_physical = None

#jg.runtime.benchmark_name = "gaussian_bumps_test_cases"

#
# Switch to Gaussian Bump since something weired goes on with PVD
# Maybe inconsistent Div field
#
#jg.runtime.benchmark_name = "gaussian_bumps_pvd"
jg.runtime.benchmark_name = "gaussian_bump"

#
# Compute error
#
jg.runtime.compute_errors = 0

jg.runtime.f_sphere = 0
jg.compile.mode = "release"

#jg.runtime.gravitation= 1
#jg.runtime.sphere_rotating_coriolis_omega = 1
#jg.runtime.h0 = 1
#jg.runtime.plane_domain_size = 1

jg.runtime.viscosity = 0.0


jg.unique_id_filter = [
        'compile',
        'parallelization',
        'runtime.galewsky_params',
        'runtime.normal_modes_params',
        'runtime.simparams',
        'runtime.timestepping_method',
        'runtime.timestepping_order',
        'runtime.semi_lagrangian',
    ]

#####################################################
#####################################################
#####################################################

ref_ts_method = 'ln_erk'
ref_ts_order = 4

# Parameters to be tested for SDC
listNodeTypes = ['RADAU-RIGHT', 'LOBATTO']
listNNodes = [ 3, 4, 5 ]
listNIters = [ 1, 2, 3 ]

ref_ts_size = 2
timestep_size_min = 64
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

jg.runtime.max_simulation_time = timestep_size_min*512

#####################################################
#####################################################
#####################################################

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

#
# Reference solution
#
jg.runtime.rexi_method = None
jg.runtime.timestepping_method = ref_ts_method
jg.runtime.timestepping_order = ref_ts_order
jg.runtime.timestepping_order2 = ref_ts_order
jg.runtime.timestep_size = ref_ts_size

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id

#
# Create job scripts
#

for nodeType in listNodeTypes:
    for nNodes in listNNodes:
        for nIter in listNIters:
            for jg.runtime.timestep_size in timestep_sizes:
                
                jg.runtime.timestepping_method = 'ln_imex_sdc'

                # Generate SDC parameter dictionary
                jg.runtime.paramsSDC = getSDCSetup(nNodes, nodeType, nIter)
                
                # Create SDC specific runtime attributes
                jg.runtime.init_phase = True
                jg.runtime.nodeType = nodeType
                jg.runtime.nNodes = nNodes
                jg.runtime.nIter = nIter
                jg.runtime.init_phase = False

                # Each sweep should give one order of accuracy (initial sweep type = copy)
                jg.runtime.timestepping_order = nIter

                if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
                    print("simtime: "+str(jg.runtime.max_simulation_time))
                    print("timestep_size: "+str(jg.runtime.timestep_size))
                    raise Exception("Invalid time step size (not remainder-less dividable)")

                jg.gen_jobscript_directory()
