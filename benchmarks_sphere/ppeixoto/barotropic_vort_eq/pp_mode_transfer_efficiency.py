#! /usr/bin/env python3

import sys

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D



from mule.plotting.Plotting import *
from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import modes_experiment as mexp

if len(sys.argv) > 1:
	experiment_file = sys.argv[1]
else:
	print("")
	print("Usage:")
	print("")
	print("	"+sys.argv[0]+" [experiment file.pckl]")
	print("")
	sys.exit(1)


#get experiment setup
filename = experiment_file #"mode_setup_1.pckl"
obj=mexp.load_file(filename)
exp_codes=obj.codes
basebanchname="barotropic_vort_modes_"
exp_codes_bnames=[basebanchname+s for s in exp_codes]
print("Benchmarks under investigation:")
print(exp_codes_bnames)

#get jobs data
j = JobsData('./job_bench_*', verbosity=1)
#print(j)
jobs = j.get_jobs_data()
#print(jobs)
main_tag = 'runtime.benchmark_name'
alpha_tag = 'output.benchmark_barotropic_vort_modes.0.ampl' #mode 0 is used as reference
jobs_flat =[]
jobs_raw =[]
alphas = []
job_dirs = []
for key, job in jobs.items():
	d = job.get_flattened_data()
	r = job.get_job_raw_data()
	if d[main_tag] in exp_codes_bnames:
		jobs_flat.append(d)
		jobs_raw.append(r)
		a = d['output.benchmark_barotropic_vort_modes.0.ampl']
		alphas.append(a)
		dir = d['jobgeneration.job_dirpath']
		job_dirs.append(dir)
		

for i in range(len(alphas)):
	print()
	print("Post-processing (alpha, dir):",	alphas[i], job_dirs[i])

	#List all parameters of job
	jd_flat = jobs_flat[i]
	#for key in jd_flat:
	#	print(key, '->', jd_flat[key])		
	jd_raw = jobs_raw[i]

	output=jd_raw['output']
	runtime=jd_raw['jobgeneration']
	runtime=runtime['runtime']
	code=output['benchmark_barotropic_vort_modes.code']
	maxmodes = int(output["benchmark_barotropic_vort_modes.maxmodes"])

	nmodes=[]
	mmodes=[]
	ampls=[]

	print("Initial conditions")
	print("Mode, n, m, amplitude")
	for i in range(maxmodes):
		mode="benchmark_barotropic_vort_modes."+str(i)+"."
		imode=i
		nmode=int(output[mode+"nmode"])
		mmode=int(output[mode+"mmode"])
		ampl=float(output[mode+"ampl"])
		print(i, nmode, mmode, ampl)
		nmodes.append(nmode)
		mmodes.append(mmode)
		ampls.append(ampl)

	evol = mexp.evol(jd_flat['runtime.p_job_dirpath'])

	nout_shell_min = 3
	nout_shell_max = 6
	evol.set_out_shells(nout_shell_min,nout_shell_max)
	
	evol.plot(code, "mode_evol.pdf")

	filename_shell = "shell_evol_n"+str(nout_shell_min)+"-"+str(nout_shell_max)+".pdf"
	evol.plot_shells(code, filename_shell)

	#exit(1)

