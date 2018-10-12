#! /usr/bin/env python3

import sys
import math

from SWEET import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

groups = ['runtime.timestepping_method']


tagnames_y = [
	'sphere_data_diff_prog_h.norm_l1',
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
			{
				'tagname_x': 'output.simulation_benchmark_timings.main_timestepping',
				'xlabel': "Wallclock time (seconds)",
				'ylabel': tagname_y,
				'title': 'Wallclock time vs. error',
				'xscale': 'log',
				'yscale': 'log',
			},
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

		d = JobsData_GroupsDataTable(
				job_groups,
				tagname_x,
				tagname_y,
				data_filter = data_filter
			)
		fileid = "output_table_"+tagname_x.replace('.', '-').replace('_', '-')+"_vs_"+tagname_y.replace('.', '-').replace('_', '-')

		print("Data table:")
		d.print()
		d.write(fileid+".csv")


		"""
		Plotting format
		"""

		d = JobsData_GroupsPlottingScattered(
				job_groups,
				tagname_x,
				tagname_y,
				data_filter = data_filter
			)

		fileid = "output_plotting_"+tagname_x.replace('.', '-').replace('_', '-')+"_vs_"+tagname_y.replace('.', '-').replace('_', '-')

		p = Plotting()
		p.plot_scattered(
				data_plotting = d.data,
				xlabel = xlabel,
				ylabel = ylabel,
				title = title,
				xscale = xscale,
				yscale = yscale,
				outfile=fileid+".pdf",
			)

		print("Data plotting:")
		d.print()
		d.write(fileid+".csv")

		print("Info:")
		print("	NaN: Errors in simulations")
		print("	None: No data available")
