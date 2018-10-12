#! /usr/bin/env python3

from SWEET import *
from SWEETPostprocessingPlotting import *
from SWEETPostprocessingJobsData import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys

if len(sys.argv) <= 1:
	print("")
	print("Usage:")
	print("	"+sys.argv[0]+" [tag name]")
	print("")
	sys.exit(1)


groups = ['runtime.timestepping_method']
tagname_y = sys.argv[1]


#params_tagname_x = ['runtime.timestep_size', 'TODO: Wallclocktime']
params_tagname_x = ['runtime.timestep_size']

for tagname_x in params_tagname_x:

	j = SWEETPostprocessingJobsData('./job_bench_*', verbosity=10)

	print("")
	print("*"*80)
	print("Groups:")
	job_groups = j.create_groups(groups)
	for g in job_groups:
		print(g)

	# Filter out errors beyond this value!
	def data_filter(x, y, jobdata):
		x = float(x)
		y = float(y)

		if y > 9.0:
			return True

		if False:
			if x > 30 and x < 100:
				return True

			if x < 20:
				return True

		return False

	data_table = j.create_data_table_float(
			groups,
			tagname_x,
			tagname_y,
			data_filter = data_filter
		)

	print("Data table:")
	j.print_data_table(data_table)

	data_plotting = j.create_data_plotting_float(
			groups,
			tagname_x,
			tagname_y,
			data_filter = data_filter
		)

	p = SWEETPostprocessingPlotting()
	p.plot_scattered(
			data_plotting=data_plotting,
			xlabel="Timestep size",
			ylabel="L1 error",
			title = "Timestep size vs. Error",
			xscale="log",
			yscale="log",
			outfile="output_dt_vs_error_l1.pdf",
		)

