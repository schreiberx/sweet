#! /usr/bin/env python3

import sys
import math

from SWEET import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

groups = ['runtime.timestepping_method', 'runtime.h']


tagnames_y = [
	'sphere_data_diff_prog_h.norm_l1',
]



j = JobsData('./job_bench_*', verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for g in job_groups:
	print(g)

# Filter out errors beyond this value!
def data_filter(group_id, x, y):
	x = float(x)
	y = float(y)

	if 'timings' in tagname_x:
		# Filter out NaNs for wallclock time studies
		# NaNs require significantly more computation time
		if math.isnan(y):
			return True
	# No filter!
	# No filter!
	# No filter!
	return False

	if y > 9.0:
		return True

	if False:
		if x > 30 and x < 100:
			return True

		if x < 20:
			return True

	return False


for tagname_y in tagnames_y:

	params = [
			{
				'tagname_x': 'runtime.timestep_size',
				'xlabel': "Timestep size",
				'ylabel': tagname_y,
				'title': 'Timestep size vs. error',
				'xscale': 'log',
				'yscale': 'log',
			},
"""
			{
				'tagname_x': 'output.simulation_benchmark_timings.main_timestepping',
				'xlabel': "Wallclock time (seconds)",
				'ylabel': tagname_y,
				'title': 'Wallclock time vs. error',
				'xscale': 'log',
				'yscale': 'log',
			},
"""
		]


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

		"""
		Table format
		"""

		data_table = c.create_data_table_float(
				groups,
				tagname_x,
				tagname_y,
				data_filter = lambda a,b,c: data_filter(a,b,c)
			)
		fileid = "output_table_"+tagname_x.replace('.', '-').replace('_', '-')+"_vs_"+tagname_y.replace('.', '-').replace('_', '-')

		print("Data table:")
		c.print_data_table(data_table)
		c.write_data_table(data_table, fileid+".csv")


		"""
		Plotting format
		"""

		data_plotting = c.create_data_plotting_float(
				groups,
				tagname_x,
				tagname_y,
				data_filter = data_filter
			)

		fileid = "output_plotting_"+tagname_x.replace('.', '-').replace('_', '-')+"_vs_"+tagname_y.replace('.', '-').replace('_', '-')

		p = Plotting()
		p.plot_scattered(
				data_plotting = data_plotting,
				xlabel = xlabel,
				ylabel = ylabel,
				title = title,
				xscale = xscale,
				yscale = yscale,
				outfile=fileid+".pdf",
			)

		print("Data plotting:")
		c.print_data_plotting(data_plotting)
		c.write_data_plotting(data_plotting, fileid+".csv")

		print("Info:")
		print("	NaN: Errors in simulations")
		print("	None: No data available")
