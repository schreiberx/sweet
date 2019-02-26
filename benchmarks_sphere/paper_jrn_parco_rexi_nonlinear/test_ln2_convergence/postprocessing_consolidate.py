#! /usr/bin/env python3

import sys
import math

from mule_local.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

groups = ['runtime.timestepping_method', 'runtime.max_simulation_time']


tagnames_y = [
	'sphere_data_diff_prog_h.norm_l1',
	'sphere_data_diff_prog_h.norm_linf',
]


j = JobsData('./job_bench_*', verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
	print(key)

# Filter out errors beyond this value!
def data_filter(x, y, jobdata):
	x = float(x)
	y = float(y)

	return False

	# Filter out NaNs for wallclock time studies
	# NaNs require significantly more computation time
	if math.isnan(y):
		return True

	if y > 10.0:
		return True

	return False


for tagname_y in tagnames_y:
	print("*"*80)
	print("Processing tagname "+tagname_y)
	print("*"*80)

	params = [
			{
				'tagname_x': 'runtime.timestep_size',
				'xlabel': "Timestep size",
				'ylabel': tagname_y,
				'title': 'Timestep size vs. error',
				'xscale': 'log',
				'yscale': 'log',
			},
		]
	"""
	params += [
			{
				'tagname_x': 'output.simulation_benchmark_timings.main_timestepping',
				'xlabel': "Wallclock time (seconds)",
				'ylabel': tagname_y,
				'title': 'Wallclock time vs. error',
				'xscale': 'log',
				'yscale': 'log',
			},
		]

	"""

	for param in params:

		tagname_x = param['tagname_x']
		xlabel = param['xlabel']
		ylabel = param['ylabel']
		title = param['title']
		xscale = param['xscale']
		yscale = param['yscale']

		print("*"*80)
		print("Processing tag "+tagname_x)
		print("*"*80)

		#if True:
		if False:
			"""
			Table format
			"""

			d = JobsData_GroupsDataTable(
					job_groups,
					tagname_x,
					tagname_y,
					data_filter = data_filter
				)
			fileid = "output_table_"+tagname_x.replace('.', '-').replace('_', '-')+"_vs_"+tagname_y.replace('.', '-').replace('_', '-')

			print("Data table:")
			d.print()
			#d.write(fileid+".csv")


		if True:
			"""
			Plotting format
			"""

			d = JobsData_GroupsPlottingScattered(
					job_groups,
					tagname_x,
					tagname_y,
					data_filter = data_filter
				)

			for key, group_data in d.get_data_float().items():
				print("Group: "+key)
				prev_value = -1.0
				conv = '-'
				for (x, y) in zip(group_data['x_values'], group_data['y_values']):
					if prev_value >= 0:
						conv = y/prev_value
					print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv))
					prev_value = y
			

		print("Info:")
		print("	NaN: Errors in simulations")
		print("	None: No data available")
