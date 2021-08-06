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
	print(alphas[i], job_dirs[i])


	#List all parameters of job
	jd_flat = jobs_flat[i]
	#for key in jd_flat:
	#	print(key, '->', jd_flat[key])		
	jd_raw = jobs_raw[i]

	output=jd_raw['output']
	runtime=jd_raw['jobgeneration']
	runtime=runtime['runtime']

print(runtime)


#modes_experiment=mexp.load_file()
#print(modes_experiment.nmodes)
exit(1)
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



#Read output_nm files
energy_file    = jd_flat['runtime.p_job_dirpath']+"/output_spec_energy_t00000000000.00000000.txt"
enstrophy_file = jd_flat['runtime.p_job_dirpath']+"/output_spec_enstrophy_t00000000000.00000000.txt"


scalesmin = {}
scalesmax = {}
timerescale=1.0
#Remove modes with null values
exclude = ['timestamp']

df_energy=pd.read_csv(energy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
print(df_energy)
df_ens=pd.read_csv(enstrophy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
print(df_ens)

time = df_energy['timestamp']
tenergy = df_energy['TotalSum']
maxenergy = df_energy['TotalSum'].max()
eps_en=maxenergy/100
scalesmin[0]=eps_en
scalesmax[0]=maxenergy
df_energy_clean = df_energy.loc[:, (df_energy > eps_en).any(axis=0)]
print(df_energy)
#df_energy.set_index('timestamp',drop=True,inplace=True)
print(df_energy_clean)

tens = df_ens['TotalSum']
maxens = df_ens['TotalSum'].max()
eps_ens=maxens/100
scalesmin[1]=eps_ens
scalesmax[1]=maxens
df_ens_clean = df_ens.loc[:, (df_ens > eps_ens).any(axis=0)]
print(df_ens)
#df_energy.set_index('timestamp',drop=True,inplace=True)
print(df_ens_clean)

##########################################################
# Plotting starts here
##########################################################

print("*"*80)
print("*"*80)
print("*"*80)


fontsize=18
figsize=(10, 10)

fig, axs = plt.subplots(2, figsize=(10,10), sharex=True)
#plt.rc('text', usetex=True)
title="SPH Mode Nonlinear Interaction\n"+code
fig.suptitle(title)


for i, ax in enumerate(axs):
	ax.set_xscale("linear")
	ax.set_yscale("log", nonpositive='clip')
	ylim=[scalesmin[i], scalesmax[i]]
	ax.set_ylim(ylim)

for ax in axs.flat:
    ax.label_outer()

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ' and m != '':
            markers.append(m)
    except TypeError:
        pass

linestyles = ['-', '--', ':', '-.']


ncol = 2	
axs[0].set(ylabel='Energy', xlabel="time")
df_energy_clean.plot(x='timestamp', ax=axs[0])
axs[0].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)


ncol=2
axs[1].set(ylabel='Enstrophy', xlabel="time")
df_ens_clean.plot(x='timestamp', ax=axs[1])
axs[1].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

fig.subplots_adjust(right=0.7)

plt.savefig(output_filename, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)

plt.close()


