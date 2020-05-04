#! /usr/bin/env python3

import sys
import math
import copy

from mule_local.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
import sys


groups = [
    	'runtime.timestepping_method',
    	'runtime.use_robert_functions',
    ]

tagname_x = 'runtime.timestep_size'
tagname_y = 'output.errors'

j = JobsData('./job_bench_*', verbosity=0)

c = JobsDataConsolidate(j)


print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
    print(key)


errors = []
for key, jobs_data in job_groups.items():
    for job_key, job_data in jobs_data.items():
    	error_line = job_data[tagname_y]
    	error_split = error_line.split("\t")
    	if len(error_split) != 4:
    		raise Exception("Inconsistent number of elements in error output")

    	# Avoid including unstable results
    	simtime = error_split[0]
    	if abs(float(simtime.replace('simtime=', '')) - job_data['runtime.max_simulation_time']) > 1e-10:
    		continue

    	error_linf_phi = float(error_split[1].replace('error_linf_phi=', ''))

    	g = 9.80616
    	gh = 29400
    	rel_error_h = error_linf_phi / gh

    	print(job_data['jobgeneration.job_dirpath']+": rel_error_h = "+str(rel_error_h))
    	errors.append(rel_error_h)


import math
import numpy

max_error = numpy.amax(numpy.abs(errors))
print("Maximum overall error: "+str(max_error))

err_threshold = 1e-9
if max_error > err_threshold or not math.isfinite(max_error):
    raise Exception("Relative overall error exceeeds threshold "+str(err_threshold))

print("Tests successful")
