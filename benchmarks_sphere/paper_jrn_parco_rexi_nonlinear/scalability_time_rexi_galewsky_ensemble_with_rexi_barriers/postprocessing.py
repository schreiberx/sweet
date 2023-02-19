#! /usr/bin/env python3

import sys

from mule import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule.plotting.Plotting import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


bar_data = 	[
				'output.simulation_benchmark_timings.main_timestepping',
				'output.simulation_benchmark_timings.rexi_timestepping',
#				'output.simulation_benchmark_timings.rexi',
#				'output.simulation_benchmark_timings.rexi_setup',
#				'output.simulation_benchmark_timings.rexi_shutdown',
				'output.simulation_benchmark_timings.rexi_timestepping_broadcast',
				'output.simulation_benchmark_timings.rexi_timestepping_solver',
				'output.simulation_benchmark_timings.rexi_timestepping_reduce',
#				'output.simulation_benchmark_timings.rexi_timestepping_miscprocessing',
		]



def arraylist_swap_cols(data, cola, colb):
	data_new = [[None for i in range(len(data[0]))] for j in range(len(data))]
	for j in range(len(data)):
		for i in range(len(data[0])):
			data_new[j][i] = data[j][i]

	for j in range(len(data)):
		tmp = data_new[j][cola]
		data_new[j][cola] = data_new[j][colb]
		data_new[j][colb] = tmp

	return data_new



def arraylist_remove_last_col(data):
	data_new = [[None for i in range(len(data[0])-1)] for j in range(len(data))]
	for j in range(len(data)):
		for i in range(len(data[0])-1):
			data_new[j][i] = data[j][i]
	return data_new


num_ensembles = 10

ensemble_data = []
for ensemble_id in range(num_ensembles):

	j = JobsData('job_bench_*ensemble'+str(ensemble_id).zfill(2)+'*', verbosity=100)

	c = JobsDataConsolidate(j)

	d = JobsData_DataTable(
			j,
			'parallelization.num_ranks',
			bar_data
		)

	data = d.get_data_float()
	if True:
	#if False:
		"""
		Add last column 'nl_timestepping'
		"""
		data_new = [[None for i in range(len(data[0])+1)] for j in range(len(data))]
		for j in range(len(data)):
			for i in range(len(data[0])):
				data_new[j][i] = data[j][i]

		keya = bar_data.index('output.simulation_benchmark_timings.main_timestepping')+1
		keyb = bar_data.index('output.simulation_benchmark_timings.rexi_timestepping')+1

		data_new[0][-1] = 'nl_timestepping'
		for j in range(1, len(data)):
			data_new[j][-1] = data[j][keya] - data[j][keyb]

		data = data_new

		data = arraylist_swap_cols(data, 1, -1)
		data = arraylist_remove_last_col(data)

	if True:
		"""
		Transpose
		"""
		data_new = [[None for i in range(len(data))] for j in range(len(data[0]))]
		for j in range(len(data)):
			for i in range(len(data[0])):
				data_new[i][j] = data[j][i]

		data = data_new
		
	ensemble_data += [data]

	for d in data:
		print(d)
	print("")


#
# Merge ensembles togehter
#
sx = len(ensemble_data[0][0])
sy = len(ensemble_data[0])
print("Dims: "+str(sx)+" x "+str(sy))

accum_data = [[None for i in range(sx)] for j in range(sy)]

# copy header
for iy in range(1, sy):
	accum_data[iy][0] = ensemble_data[0][iy][0]

for ix in range(1, sx):
	accum_data[0][ix] = ensemble_data[0][0][ix]
	accum_data[0][ix] = str(accum_data[0][ix]).replace('output.simulation_benchmark_timings.', '')

for iy in range(1, sy):
	for ix in range(1, sx):
		accum_data[iy][ix] = ensemble_data[0][iy][ix]
	accum_data[iy][0] = str(accum_data[iy][0]).replace('output.simulation_benchmark_timings.', '')


ensemble_range = range(1, len(ensemble_data))
for ie in ensemble_range:
	e = ensemble_data[ie]
	for ix in range(1, sx):
		for iy in range(1, sy):
			accum_data[iy][ix] = min(accum_data[iy][ix], e[iy][ix])



for ix in range(1, sx):
	accum_data[0][ix] = accum_data[0][ix]+" MPI ranks"

p = Plotting_Bars()

p.gen_plot_from_tabledata(
		accum_data,
		#xlabel="Number of ranks",
		ylabel="Wallclock time",
		yscale="log",
		title="Breakdown of NL-REXI wallclock times (min filter over "+str(num_ensembles)+" ensembles)",
		#subtitle=str(num_ensembles)+" ensembles using min. filter",
		outfile="output_breakdown_wallclocktimes.pdf",
		annotate_bars_with_values=True,
		annotate_bars_with_labels=True,
		legend=False,
		ylim=(1e-12, 200),
		filled_bars=False,
	)
