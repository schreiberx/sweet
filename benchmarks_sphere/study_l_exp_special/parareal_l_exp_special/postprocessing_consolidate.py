#! /usr/bin/env python3

import sys
import math
##import shutil
##import pickle

from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule_local.postprocessing.PInT_Errors import *

from matplotlib.ticker import MaxNLocator

sys.path.append('../')
import pretty_plotting as pp
sys.path.pop()

mule_plotting_usetex(False)


#### TODO
##### add command-line argument to choose between
##### - plot vs iteration
##### - plot fixed iteration vs parareal_coarse_timestep_size

"""
Postprocessing script
Plot parareal errors and residuals along iterations for simulations with given parameter values
plot_type : "iteration" or "timestep"

Run:
        ./postprocessing_consolidate plot_type var1 "name_var1" "vals_var1" var2 "names_var2" "vals_var2"
Example
        ./postprocessing_consolidate plot_type var1 runtime.parareal_coarse_slices 5 var2 runtime.parareal_coarse_timestep_size 120 240 360
"""


groups = [
             ####'runtime.timestepping_method',
             'runtime.parareal_coarse_timestepping_method',
             'runtime.parareal_coarse_timestepping_order',
             'runtime.parareal_coarse_slices'
         ]



## get filter params
filter_params = {}

for i in range(1, len(sys.argv)):
    if sys.argv[i][:3] == "var":
        var = sys.argv[i + 1]
        j = i + 2
        vals = []
        while True:
            if j == len(sys.argv):
                break
            if sys.argv[j][:3] == "var":
                break
            vals.append(sys.argv[j])
            j += 1
        filter_params[var] = vals


plot_type = sys.argv[1]
error_type = sys.argv[2]
plot_timings = int(sys.argv[3])
if not (plot_type == "iteration" or plot_type == "timestep"):
    raise Exception("Invalid plot_type")
if plot_type == "iteration" and "runtime.iteration" in filter_params.keys():
    raise Exception("'runtime.iteration' should not be set as a filter parameter")
if plot_type == "timestep" and "runtime.iteration" not in filter_params.keys():
    raise Exception("'runtime.iteration' should be set as a filter parameter")

if plot_type == "iteration":
    groups.append('runtime.parareal_coarse_timestep_size')
if plot_type == "timestep":
    groups.append('runtime.iteration')

if not (error_type == "physical" or error_type == "spectral"):
    raise Exception("Wrong error_type: " + error_type);

print(filter_params)

if error_type == "physical":
    tagnames_y = [
            'parareal_errors.prog_phi_pert.t1.0.norm_l1',
            'parareal_errors.prog_phi_pert.t1.0.norm_l2',
            'parareal_errors.prog_phi_pert.t1.0.norm_linf',
            'parareal_errors.prog_vrt.t1.0.norm_l1',
            'parareal_errors.prog_vrt.t1.0.norm_l2',
            'parareal_errors.prog_vrt.t1.0.norm_linf',
            'parareal_errors.prog_div.t1.0.norm_l1',
            'parareal_errors.prog_div.t1.0.norm_l2',
            'parareal_errors.prog_div.t1.0.norm_linf'
    ]
elif error_type == "spectral":
    tagnames_y = [
            'parareal_errors.spec.prog_phi_pert.t1.0.norm_linf_rnorm16',
            'parareal_errors.spec.prog_phi_pert.t1.0.norm_linf_rnorm32',
            'parareal_errors.spec.prog_vrt.t1.0.norm_linf_rnorm16',
            'parareal_errors.spec.prog_vrt.t1.0.norm_linf_rnorm32',
            'parareal_errors.spec.prog_div.t1.0.norm_linf_rnorm16',
            'parareal_errors.spec.prog_div.t1.0.norm_linf_rnorm32'
    ]

####filter_params = {
####                     'runtime.parareal_coarse_slices' : [5]
####                }


j = JobsData('./job_bench_*', verbosity=0)
c = JobsDataConsolidate(j)
c_pint = JobsDataPInTConsolidate(j)



## Small hack for plotting error vs iteration using JobsData_GroupsPlottingScattered
## and also for plotting error at fixed iteration (in function e.g. of coarse timestep size)
## For each job job_name, make one copy job_name_iterXXX per iteration
## modifying its jobgeneration.pickle:
##   - Modify field jobgeneration.job_unique_id
##   - Modify field jobgeneration.job_dirpath
##   - Create field runtime.iteration
## Then rerun j = JobsData(...)

## Unique job per iteration
c_pint.create_copy_for_each_iteration();

## read jobs again
j = JobsData('./job_bench_*_iter*', verbosity=0)
c = JobsDataConsolidate(j)
jobs = j.get_flattened_data()

## keep simulations with given parameters
for var, val in filter_params.items():
	jobs_copy = jobs.copy();
	for job in jobs_copy.keys():
		if str(jobs[job][var]) not in val:
			del(jobs[job])

print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
	print(key)

pint_errors = PInT_Errors(
                           job_directories = None,
                           pint_type = "parareal",
                           geometry = "sphere",
                           precomputed_errors = True,
                           file_type = "csv",
                           ref_type = "ref",
                           error_type = error_type,
                           jobdir_pattern = None
              )


err = pint_errors.get_pint_errors()

for tagname_y in tagnames_y:

	err_title = "";
	if error_type == "physical":
		err_title = " (physical)"
	elif error_type == "spectral":
		## find rnorm value
		rnorm = tagname_y[tagname_y.find("rnorm") + 5:]
		print("Rnorm", rnorm);
		err_title = " (spectral, " + r'$R_{norm} = $' + r'${}$'.format(rnorm) + ")"

	if plot_type == "iteration":
		params = []
		params += [
				{
					'tagname_x': 'runtime.iteration',
					'xlabel': "Iteration",
					'ylabel': pp.latex_pretty_names[tagname_y],
					'title': 'Iteration vs. error' + err_title,
					'xscale': 'linear',
					'yscale': 'log',
				},
			]

		if plot_timings:
			raise Exception("Plot timings not implemented.")
			params += [
					{
						'tagname_x': 'output.simulation_benchmark_timings.main_timestepping',
						'xlabel': "Wallclock time (seconds)",
						'ylabel': pp.latex_pretty_names[tagname_y],
						'title': 'Wallclock time vs. error' + err_title,
						'xscale': 'log',
						'yscale': 'log',
					},
				]
	elif plot_type == "timestep":
		params = []
		params += [
				{
					'tagname_x': 'runtime.parareal_coarse_timestep_size',
					'xlabel': "Coarse timestep size",
					'ylabel': pp.latex_pretty_names[tagname_y],
					'title': 'Coarse timestep size vs. error' + err_title,
					'xscale': 'log',
					'yscale': 'log',
				},
			]

		if plot_timings:
			raise Exception("Plot timings not implemented.")
			params += [
					{
						'tagname_x': 'output.simulation_benchmark_timings.main_timestepping',
						'xlabel': "Wallclock time (seconds)",
						'ylabel': pp.latex_pretty_names[tagname_y],
						'title': 'Wallclock time vs. error' + err_title,
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

				if 'prog_h' in tagname_y:
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
			for var, val in filter_params.items():
				fileid += "_" + pp.get_abbrev_param_name(var);
				for v in val:
					fileid += '_' + str(v)
				fileid += "_"
			## remove final "_"
			if filter_params:
				fileid = fileid[:-1]

			if True:
				#
				# Proper naming and sorting of each label
				#

				# new data dictionary
				data_new = {}
				for key, data in d.data.items():



					# generate nice tex label
					#data['label'] = pp.get_pretty_name(key)
					data['label'] = key #pp.get_pretty_name(key)

					key_new = pp.get_pretty_name_order(key)+'_'+key

					# copy data
					data_new[key_new] = copy.copy(data)
 

				# Copy back new data table
				d.data = data_new

			p = Plotting_ScatteredData()


			def fun(p):
				from matplotlib import ticker
				from matplotlib.ticker import FormatStrFormatter

				if plot_type == "iteration":
					p.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
				elif plot_type == "timestep":
					###pass
					plt.tick_params(axis='x', which='minor')
					p.ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
					p.ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))

					p.ax.xaxis.set_minor_locator(ticker.LogLocator(subs=[1.5, 2.0, 3.0, 5.0]))


				for tick in p.ax.xaxis.get_minor_ticks():
					tick.label.set_fontsize(8) 


				plt.tick_params(axis='y', which='minor')
				p.ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1e"))
				p.ax.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))
 
				###p.ax.yaxis.set_minor_locator(ticker.LogLocator(subs=[1.5, 2.0, 3.0, 5.0]))

				for tick in p.ax.yaxis.get_minor_ticks():
					tick.label.set_fontsize(6) 



			annotate_text_template = "{:.1f} / {:.3f}"
			p.plot(
					data_plotting = d.get_data_float(),
					xlabel = xlabel,
					ylabel = ylabel,
					title = title,
					xscale = xscale,
					yscale = yscale,
					#annotate = True,
					#annotate_each_nth_value = 3,
					#annotate_fontsize = 6,
					#annotate_text_template = annotate_text_template,
					legend_fontsize = 8,
					grid = True,
					outfile = fileid+".pdf",
					lambda_fun = fun,
				)

			print("Data plotting:")
			d.print()
			d.write(fileid+".csv")

		print("Info:")
		print("	NaN: Errors in simulations")
		print("	None: No data available")
