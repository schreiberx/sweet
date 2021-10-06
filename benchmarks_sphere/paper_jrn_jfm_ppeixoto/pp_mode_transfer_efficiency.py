#! /usr/bin/env python3

import sys
import os

import matplotlib
#matplotlib.use('Agg')
matplotlib.use('TkAgg')

matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import matplotlib.ticker as mtick

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

if len(sys.argv) > 2:
	print(sys.argv[2])
	if int(sys.argv[2])>0:
		plots = True
	else:
		plots = False
else:
	plots = True


#get experiment setup
print("Input experiment file:", experiment_pckl_file)
pckl_file = experiment_pckl_file #"mode_setup_1.pckl"
pckl_filename = os.path.basename(pckl_file)
basedir = os.path.dirname(pckl_file)

exp_tag = os.path.splitext(pckl_filename)[0]

print("Experiment Dir:", basedir)
print("Experiment Tag:", exp_tag )
print()

#Load experiment info
obj = mexp.load_file(pckl_file)
exp_codes = obj.codes
basebanchname = "barotropic_vort_modes_"
exp_codes_bnames = [basebanchname+s for s in exp_codes]
#print("Benchmarks under investigation:")
#print(exp_codes_bnames)

title_in = "Init Modes:"
n = len(obj.nmodes)
for i in range(n):
	title_in = title_in + " ("+str(obj.nmodes[i])+";"+str(obj.mmodes[i])+")"

print(title_in)

#get jobs data
print("Extracting jobs")
j = JobsData(basedir+'/job_bench_*', verbosity=1)
jobs = j.get_jobs_data()
print("...done")

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

#print("alpha   umax  vmax" )
for key, job in sorted(jobs.items()):
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
		umax.append(float(u))
		v = d['output.benchmark_barotropic_vort_modes.vmax']
		vmax.append(float(v))
		#print(a, u, v)

df_reduced = pd.DataFrame(list(zip(alphas, umax, vmax)), columns =['Alpha', 'u_max', 'v_max'])
df_reduced = df_reduced.set_index('Alpha')
df_reduced = df_reduced.sort_index()
print(df_reduced)

#Plot velocities
if plots:
	plt.figure(figsize=(10,6), tight_layout=True)
	plt.plot(df_reduced, '-', linewidth=2)

	plt.xlabel(r" $\alpha$")
	plt.ylabel(r" Max Velocity $m/s$")
	title = title_in
	plt.title(title)
	plt.legend(title_fontsize = 13, labels=['zonal (u)' , 'meridional (v)'])

	filename_final = basedir+"/"+exp_tag+"_vel.pdf"
	print("Vel file:", filename_final)
	plt.savefig(filename_final, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)
	plt.close()

#Calculate efficiency
max_exchange_out_energy = []
max_exchange_noninit_energy = []
max_exchange_out_ens = []
max_exchange_noninit_ens = []

dominant_period = []

#output shells
print(basedir)
lim_inf=0
lim_sup=9999999
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
	#n_out_list = [5]
	#m_out_list = [3]
	nout_shell_min = 7
	nout_shell_max = 7
	out_type = "mode"
	trunc = True
	trunc_alpha = 30
	#fourier period truncation
	lim_inf=12
	lim_sup=120

elif "TC2_in5-4_7-3" in basedir:
	n_out_list = [9]
	m_out_list = [2]
	nout_shell_min = 7
	nout_shell_max = 7
	out_type = "mode"
elif "TC3_in7-3_3-1" in basedir:
	n_out_list = [7, 3]
	m_out_list = [4, 0]
	nout_shell_min = 7
	nout_shell_max = 7
	out_type = "mode"
else:
	nout_shell_min = 7
	nout_shell_max = 7
	print("Dont know this test case")
	exit(1)

#Nice Title 
title_out = "Out Mode:"
for i in range(len(n_out_list)):
	#print("(", str(arr[i*3+0]), ";", str(arr[i*3+1]), ")")
	title_out = title_out + " ("+str(n_out_list[i])+";"+str(m_out_list[i])+")"


#Loop over all experiments
for i in range(len(alphas)):
	
	if trunc:
		if alphas[i] > trunc_alpha:
			max_exchange_out_energy.append(0)
			max_exchange_noninit_energy.append(0)
			max_exchange_out_ens.append(0)
			max_exchange_noninit_ens.append(0)
			dominant_period.append(0)
			continue

	#if alphas[i] != 20:
	#	continue

	print()
	print("Post-processing (alpha, dir, umax, vmax):\n   ",	alphas[i], job_dirs[i], umax[i], vmax[i])
	
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

	max_exchange_out_energy.append(evol.max_exchange_out_energy)
	max_exchange_noninit_energy.append(evol.max_exchange_noninit_energy)
	max_exchange_out_ens.append(evol.max_exchange_out_ens)
	max_exchange_noninit_ens.append(evol.max_exchange_noninit_ens)
	
	title = title_in + " alpha = " + str(alphas[i])
	if plots:	
		#evol.plot(title, "mode_evol.pdf")
		evol.plot(title, "mode_evol.png")

	title = title + " , "+title_out
	if plots:	
		evol.plot_out(title, filename_out+"a"+str(alphas[i])+".png")

	#fourier modes
	title = " alpha = " + str(alphas[i])+" , "+title_out
	spec_enegy = evol.fourier_modes(title, filename_out+"a"+str(alphas[i])+"_spec.png", do_plot=plots, lim_inf=lim_inf, lim_sup=lim_sup)
	dominant_period.append(spec_enegy)

	#Phase analysis
	#dif_mean = evol.phase()
	#print(dif_mean['(3;0)'])

#Sanity check!
if len(alphas) > len(max_exchange_out_energy):
	print("lengths not matching!")
	exit()

#Dataframe for efficiency plot
df = pd.DataFrame(list(zip(alphas, max_exchange_out_energy, max_exchange_out_ens, dominant_period, umax, vmax)), columns =['Alpha', 'Exch_energy', 'Exch_ens', 'SpecEnerg', "u_max", "v_max"])
df = df.set_index('Alpha')
df = df.sort_index()

#truncate
if trunc:
	df = df.truncate(after=trunc_alpha)

#put in %
df = df[['Exch_energy', 'Exch_ens', 'SpecEnerg']]
df['Exch_energy']=df['Exch_energy']*100
df['Exch_ens']=df['Exch_ens']*100
print(df)

fig, ax = plt.subplots(figsize=(5,5)) #, tight_layout=True)
plot_energy = True
if plot_energy:
	ax.plot(df.index, df['Exch_energy'],  linewidth=1, color = 'blue', linestyle = '-', label='Energy') #, '--', ':'])
	ax.set_ylabel(r" $\epsilon(\alpha)$", color='blue', fontsize=16)
else:
	ax.plot(df.index, df['Exch_ens'],  linewidth=1, color = 'green', linestyle = '-', label='Enstrophy') #, '--', ':'])
	ax.set_ylabel(r" $\epsilon(\alpha)$", color='green', fontsize=16)
ax.set_xlabel(r" $\alpha$", fontsize=16)

ax.yaxis.set_major_formatter(mtick.PercentFormatter())

if lim_sup<200:
	ax2=ax.twinx()
	spec_string= '['+str(lim_inf)+','+str(lim_sup)+'] dys'
	ax2.plot(df.index, df['SpecEnerg'],  linewidth=1, color = 'red', linestyle = ':', label='Low Freq Spec\n'+spec_string) #, '--', ':'])
	ax2.set_ylabel(r'$\sqrt{PSD}$', color='red')

title = title_in+" "+title_out
plt.title(title, fontsize=14)

ax.legend(loc=2, fontsize=12)
ax2.legend(loc=4, fontsize=11)

if plot_energy:
	filename_final = basedir+"/"+exp_tag+"_"+filename_out+"_energy.png"
else:
	filename_final = basedir+"/"+exp_tag+"_"+filename_out+"_ens.png"
print("Output file:", filename_final)
plt.tight_layout()
plt.savefig(filename_final, transparent=True,  dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)
plt.close()
