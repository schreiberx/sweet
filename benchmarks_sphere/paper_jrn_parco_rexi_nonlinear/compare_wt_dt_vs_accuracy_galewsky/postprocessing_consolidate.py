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


tagnames_and_prettynames_y = {}
tagnames_and_prettynames_y['sphere_data_diff_prog_h.norm_l1'] = "$L_1$ norm surface height"
tagnames_and_prettynames_y['sphere_data_diff_prog_h.norm_l2'] = "$L_2$ norm surface height"
tagnames_and_prettynames_y['sphere_data_diff_prog_h.norm_linf'] = "$L_{\infty}$ norm surface height"

tagnames_and_prettynames_y['sphere_data_diff_prog_vort.norm_l1'] = "$L_1$ norm surface height"
tagnames_and_prettynames_y['sphere_data_diff_prog_vort.norm_l2'] = "$L_2$ norm surface height"
tagnames_and_prettynames_y['sphere_data_diff_prog_vort.norm_linf'] = "$L_{\infty}$ norm surface height"

tagnames_and_prettynames_y['sphere_data_diff_prog_div.norm_l1'] = "normalized $L_1$, divergence field"
tagnames_and_prettynames_y['sphere_data_diff_prog_div.norm_l2'] = "normalized $L_2$, divergence field"
tagnames_and_prettynames_y['sphere_data_diff_prog_div.norm_linf'] = "$L_{\infty}$ norm, divergence field"

tagnames_y = tagnames_and_prettynames_y.keys()



j = JobsData('./job_bench_*', verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
	print(key)

for tagname_y in tagnames_y:

	params = []
	params += [
			{
				'tagname_x': 'runtime.timestep_size',
				'xlabel': "Timestep size (seconds)",
				'ylabel': tagnames_and_prettynames_y[tagname_y],
				'title': 'Timestep size vs. error',
				'xscale': 'log',
				'yscale': 'log',
			},
		]
	params += [
			{
				'tagname_x': 'output.simulation_benchmark_timings.main_timestepping',
				'xlabel': "Wallclock time (seconds)",
				'ylabel': tagnames_and_prettynames_y[tagname_y],
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



		if True:
			"""
			Plotting format
			"""

			# Filter out errors beyond this value!
			def data_filter(x, y, jobdata):
				if y == None:
					return True

				x = float(x)
				y = float(y)

				if math.isnan(y):
					return True

				if 'l1' in tagname_y:
					if y > 10.0:
						return True
				elif 'l2' in tagname_y:
					if y > 0.1:
						return True
				elif 'linf' in tagname_y:
					if y > 100.0:
						return True
				else:
					raise Exception("Unknown y tag "+tagname_y)

				return False



			d = JobsData_GroupsPlottingScattered(
					job_groups,
					tagname_x,
					tagname_y,
					data_filter = data_filter
				)

			fileid = "output_plotting_"+tagname_x.replace('.', '-').replace('_', '-')+"_vs_"+tagname_y.replace('.', '-').replace('_', '-')


			if True:
				#
				# Proper naming and sorting of each label
				#

				# new data dictionary
				data_new = {}
				for key, data in d.data.items():
					# generate nice tex label
					pretty_names = {
						'ln_erk' : "$L(U)+N(U)$: RK",

						'l_irk_n_erk_ver0' : "$L(U)$: CN, $N(U)$: RK, SSv0",
						'l_irk_n_erk_ver1' : "$L(U)$: CN, $N(U)$: RK, SSv1",

						'lg_irk_lc_n_erk_ver0' : "$L_g(U)$: CN, $L_c(U)+N(U)$: RK, SSv0",
						'lg_irk_lc_n_erk_ver1' : "$L_g(U)$: CN, $L_c(U)+N(U)$: RK, SSv1",

						'l_rexi_n_erk_ver0' : "$L(U)$: REXI, $N(U)$: RK, SSv0",
						'l_rexi_n_erk_ver1' : "$L(U)$: REXI, $N(U)$: RK, SSv1",
						'lg_rexi_lc_n_erk_ver0' : "$L_g(U)$: REXI, $L_c(U)+N(U)$: RK, SSv0",
						'lg_rexi_lc_n_erk_ver1' : "$L_g(U)$: REXI, $L_c(U)+N(U)$: RK, SSv1",
						'l_rexi_n_etdrk' : "$L(U)$, $N(U)$: ETD2RK",
						'lg_rexi_lc_n_etdrk' : "$L_g(U)$, $L_c(U)+N(U)$: ETD2RK",
					}

					pretty_names_order = [
						'ln_erk',

						'l_irk_n_erk_ver0',
						'l_irk_n_erk_ver1',

						'lg_irk_lc_n_erk_ver0',
						'lg_irk_lc_n_erk_ver1',

						'l_rexi_n_erk_ver0',
						'l_rexi_n_erk_ver1',
						'lg_rexi_lc_n_erk_ver0',
						'lg_rexi_lc_n_erk_ver1',
						'l_rexi_n_etdrk',
						'lg_rexi_lc_n_etdrk',
					]
					if key not in pretty_names:
						data['label'] = key.replace('_', '\\_')
					else:
						data['label'] = pretty_names[key]

					# reduce key to it's main information
					key_new = str(pretty_names_order.index(key)).zfill(3)+'_'+key

					# copy data
					data_new[key_new] = copy.copy(data)

				# Copy back new data table
				d.data = data_new

			p = Plotting_ScatteredData()
			p.plot(
					data_plotting = d.get_data_float(),
					xlabel = xlabel,
					ylabel = ylabel,
					title = title,
					xscale = xscale,
					yscale = yscale,
					legend_fontsize = 8,
					outfile = fileid+".pdf",
				)

			print("Data plotting:")
			d.print()
			d.write(fileid+".csv")

		print("Info:")
		print("	NaN: Errors in simulations")
		print("	None: No data available")
