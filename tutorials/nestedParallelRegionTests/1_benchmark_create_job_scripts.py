#! /usr/bin/env python3


from mule import *
from mule.JobMule import *
j = JobGeneration()

j.compile.program = 'tutorialNestedParallelRegionTests'
j.unique_id_filter = ['runtime', 'compile']
j.unique_id_filter +=  ['simparams', 'benchmark', 'timestep_order', 'timestep_size', 'disc_space', 'compile']


# Request dedicated compile script
j.compilecommand_in_jobscript = False

j.parallelization.max_wallclock_seconds = "00:01:00"


for i in range(1, j.platform_resources.num_cores_per_socket+1):

	pspace = JobParallelizationDimOptions()
	pspace.num_cores_per_rank = j.platform_resources.num_cores_per_socket
	pspace.num_threads_per_rank = i
	pspace.num_ranks = 1

	j.setup_parallelization([pspace])
	j.parallelization.core_oversubscription = False

	j.parallelization.core_affinity = 'compact'
	j.gen_jobscript_directory()

#	j.parallelization.core_affinity = 'scatter'
#	j.gen_jobscript_directory('job_'+j.getUniqueID())


j.write_compilecommands()
