#! /usr/bin/env python3

import sys
import os

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
#import seaborn as sns

#sns.set_style('darkgrid') # darkgrid, white grid, dark, white and ticks
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
plt.rc('legend', fontsize=13)    # legend fontsize
plt.rc('font', size=13)          # controls default text sizes

import numpy as np
import pandas as pd



from mule.plotting.Plotting import *
from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import modes_experiment as mexp

if len(sys.argv) > 1:
	experiment_pckl_file = sys.argv[1]
else:
	print("")
	print("Usage:")
	print("")
	print("	"+sys.argv[0]+" [experiment file.pckl]")
	print("")
	sys.exit(1)

#Sweet!
sweetroot = os.environ.get('MULE_SOFTWARE_ROOT')

#get experiment setup
print("Input experiment file:", experiment_pckl_file)
pckl_file = experiment_pckl_file #"mode_setup_1.pckl"
pckl_filename = os.path.basename(pckl_file)
basedir = os.path.dirname(pckl_file)

exp_tag = os.path.splitext(pckl_filename)[0]

print("Experiment Dir:", basedir)
print("Experiment tag:", exp_tag )
print()

#Load experiment info
obj = mexp.load_file(pckl_file)
exp_codes = obj.codes
basebanchname = "barotropic_vort_modes_"
exp_codes_bnames = [basebanchname+s for s in exp_codes]
print("Benchmarks under investigation:")
print(exp_codes_bnames)



#get jobs data
j = JobsData(basedir+'/job_bench_*', verbosity=1)
print(j)
print("Extracting jobs")
jobs = j.get_jobs_data()

#print(jobs)
main_tag = 'runtime.benchmark_name'
alpha_tag = 'output.benchmark_barotropic_vort_modes.0.ampl' #mode 0 is used as reference
jobs_flat =[]
jobs_raw =[]
alphas = []
umax = []
vmax = []
job_dirs = []
#print("jobs extracted")
#print(jobs)

print("alpha   umax  vmax" )
for key, job in jobs.items():
	d = job.get_flattened_data()
	r = job.get_job_raw_data()

	outfile_path = d["jobgeneration.p_job_stdout_filepath"]
	
	#In case this was generated in a server, fix the path
	outfile_path = outfile_path.replace(d["jobgeneration.sweetroot"], sweetroot)
	
	outfile_exists = os.path.isfile(outfile_path)	
	if d[main_tag] in exp_codes_bnames and outfile_exists:
		jobs_flat.append(d)
		jobs_raw.append(r)
		a = d['output.benchmark_barotropic_vort_modes.0.ampl']
		alphas.append(float(a))
		dir = d['jobgeneration.job_dirpath']
		job_dirs.append(dir)
		u = d['output.benchmark_barotropic_vort_modes.umax']
		umax.append(u)
		v = d['output.benchmark_barotropic_vort_modes.vmax']
		vmax.append(v)
		print(a, u, v)


max_exchange_out_energy = []
max_exchange_noninit_energy = []
max_exchange_out_ens = []
max_exchange_noninit_ens = []

dominant_period = []

#output shells
print(basedir)

if "TC1_in2_3_4" in basedir:
	nout_shell_min = 5
	nout_shell_max = 6
	n_out_list = [5]
	m_out_list = [1]
	out_type = "shell"
elif "TC2_in3-3_4-3_2-1" in basedir:
	nout_shell_min = 5
	nout_shell_max = 6
	n_out_list = [5]
	m_out_list = [1]
	out_type = "mode"
elif "TC2_in5-4_3-1_7-3" in basedir:
	n_out_list = [9]
	m_out_list = [2]
	nout_shell_min = 7
	nout_shell_max = 7
	out_type = "mode"
elif "TC3xx" in basedir:
	nout_shell_min = 7
	nout_shell_max = 7
	n_out_list = [7,5]
	m_out_list = [3,4]
	out_type = "mode"
elif "TC4xx" in basedir:
	nout_shell_min = 7
	nout_shell_max = 7
	n_out_list = [7,5]
	m_out_list = [3,4]
	out_type = "mode"
elif "TC5xx" in basedir:
	nout_shell_min = 5
	nout_shell_max = 7
	n_out_list = [7,9, 5]
	m_out_list = [3,2, 5]
	out_type = "mode"
else:
	nout_shell_min = 7
	nout_shell_max = 7
	print("Dont know this test case")
	exit(1)

for i in range(len(alphas)):
	i=70
#for i in range(4):
	print()
	print("Post-processing (alpha, dir, umax, vmax):",	alphas[i], job_dirs[i], umax[i], vmax[i])
	
	#List all parameters of job
	jd_flat = jobs_flat[i]
	#for key in jd_flat:
	#	print(key, '->', jd_flat[key])		
	jd_raw = jobs_raw[i]

	output = jd_raw['output']
	jobgen = jd_raw['jobgeneration']
	runtime = jobgen['runtime']

	job_sweet_dir = jobgen['jobgeneration']['sweetroot']
	code = output['benchmark_barotropic_vort_modes.code']
	dir_path = runtime['p_job_dirpath']

	#fix path, in case generated in server
	dir_path = dir_path.replace(job_sweet_dir, sweetroot)

	evol = mexp.evol(dir_path)
	
	if out_type == "shell":
		filename_out = evol.set_out_shells(nout_shell_min, nout_shell_max) 
	elif out_type == "mode":
		filename_out = evol.set_out_modes(n_out_list, m_out_list) 
	else:
		print("Error: please set correct out type")
		exit(1)

	#fourier modes
	evol.fourier_modes()
	dominant_period.append(evol.large_periods_energy[0])

	max_exchange_out_energy.append(evol.max_exchange_out_energy)
	max_exchange_noninit_energy.append(evol.max_exchange_noninit_energy)
	max_exchange_out_ens.append(evol.max_exchange_out_ens)
	max_exchange_noninit_ens.append(evol.max_exchange_noninit_ens)
	
	
	#print(evol.df_energy_clean)
	#print(evol.df_ens_clean)
	arr = code.split('_')
	title = "Init Modes:"
	n = int(arr[0])
	arr = arr[1:]
	for i in range(n):
		#print("(", str(arr[i*3+0]), ";", str(arr[i*3+1]), ")")
		title = title + " ("+str(arr[i*3+0])+";"+str(arr[i*3+1])+")"

	evol.plot(title, "mode_evol.pdf")
	
	title = title + "\n Out Modes:"
	for i in range(len(n_out_list)):
		#print("(", str(arr[i*3+0]), ";", str(arr[i*3+1]), ")")
		title = title + " ("+str(n_out_list[i])+";"+str(m_out_list[i])+")"

	evol.plot_out(title, filename_out+ ".pdf")
	exit()

#print(alphas, max_exchange_out_energy, max_exchange_out_ens)
df = pd.DataFrame(list(zip(alphas, max_exchange_out_energy, max_exchange_out_ens, dominant_period, umax, vmax)), columns =['Alpha', 'Exch_energy', 'Exch_ens', "Period", "u_max", "v_max"])
df = df.set_index('Alpha')

df = df.sort_index()
print(df)

out_energy = filename_out # "out_n"+str(nout_shell_min)+"_"+str(nout_shell_max)

#df_periods = df[['Period']]
#df_periods.plot( title=exp_tag+" "+out_energy)

df = df[['Exch_energy', 'Exch_ens']]
#df['Exch_energy']=df['Exch_energy']/df.index
#df['Exch_ens']=df['Exch_ens']/df.index

plt.figure(figsize=(10,6), tight_layout=True)
#plotting
plt.plot(df, '-', linewidth=2)
#customization
#plt.xticks([2017, 2018, 2019, 2020, 2021])
plt.xlabel(r"Alpha ($\alpha$)")
plt.ylabel(r"Efficiency")
#plt.ylabel(r"Efficiency/$\alpha$")
plt.title(exp_tag+" "+out_energy)
plt.legend(title_fontsize = 13, labels=['Energy' , 'Enstrophy'])

#df.plot( title=exp_tag+" "+out_energy)
#plt.show()

filename_final = basedir+"/"+exp_tag+"_"+out_energy+".pdf"
print("Output file:", filename_final)
#plt.show()
plt.savefig(filename_final, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)

plt.close()

plt.show()