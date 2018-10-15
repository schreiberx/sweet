#! /usr/bin/env python3

import sys
import math
import copy

from SWEET import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

sys.path.append('../')
import pretty_plotting as pp
sys.path.pop()

groups = [
		'runtime.timestepping_method',
		'runtime.h',
	]


tagnames_y = [
	'sphere_data_diff_prog_h.res_norm_l1',
	'sphere_data_diff_prog_h.res_norm_l2',
	'sphere_data_diff_prog_h.res_norm_linf',
]



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
				'ylabel': pp.latex_pretty_names[tagname_y],
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
				'ylabel': pp.latex_pretty_names[tagname_y],
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
			d.write(fileid+".csv")


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

				# inconsistency of REXI
				# finding right coefficients is part of future work
				if x < 20:
					return True

				if 'prog_h' in tagname_y:
					return False

					if 'l1' in tagname_y:
						if y > 1e1:
							print("Sorting out L1 data "+str(y))
							return True
					elif 'l2' in tagname_y:
						if y > 1e1:
							print("Sorting out L2 data "+str(y))
							return True
					elif 'linf' in tagname_y:
						if y > 1e2:
							print("Sorting out Linf data "+str(y))
							return True
					else:
						raise Exception("Unknown y tag "+tagname_y)

				else:
					print("TODO")

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
					# reduce key to it's main information
					key_new = key
					pos = key_new.find('.0')
					if pos >= 0:
						key_new = key_new[0:pos]
					key_new = key_new.replace('lg_rexi_lc_n_erk_ver0_', '')
					key_new = key_new.replace('lg_rexi_lc_n_erk_ver1_', '')
					key_new = key_new.replace('lg_irk_lc_n_erk_ver0_', '')
					key_new = key_new.replace('lg_irk_lc_n_erk_ver1_', '')

					# generate nice tex label
					label_new = pp.get_pretty_name(key_new)
					label_new = '$h_0 = '+label_new+' m$'
					data['label'] = label_new

					# Finally, care about proper key name for correct sorting
					key_new = key_new.zfill(8)

					# copy data
					data_new[key_new] = copy.copy(data)

				# Copy back new data table
				d.data = data_new

			p = Plotting_ScatteredData()


			def fun(p):
				from matplotlib import ticker
				from matplotlib.ticker import FormatStrFormatter

				plt.tick_params(axis='x', which='minor')
				p.ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
				p.ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))

				p.ax.xaxis.set_minor_locator(ticker.LogLocator(subs=[1.5, 2.0, 3.0, 5.0]))

				for tick in p.ax.xaxis.get_minor_ticks():
					tick.label.set_fontsize(8) 


				plt.tick_params(axis='y', which='minor')
				p.ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1e"))
				p.ax.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))
 
				p.ax.yaxis.set_minor_locator(ticker.LogLocator(subs=[1.5, 2.0, 3.0, 5.0]))

				for tick in p.ax.yaxis.get_minor_ticks():
					tick.label.set_fontsize(6) 



			p.plot(
					data_plotting = d.get_data_float(),
					xlabel = xlabel,
					ylabel = ylabel,
					title = title,
					xscale = xscale,
					yscale = yscale,
					outfile = fileid+".pdf",
					lambda_fun = fun,
					xlim = (10, 1000),
				)

			print("Data plotting:")
			d.print()
			d.write(fileid+".csv")

		print("Info:")
		print("	NaN: Errors in simulations")
		print("	None: No data available")
