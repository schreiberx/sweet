#! /usr/bin/env python3

import sys

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D


from mule.postprocessing.JobData import *

if len(sys.argv) > 2:
	output_filename = sys.argv[1]
else:
	print("")
	print("Usage:")
	print("")
	print("	"+sys.argv[0]+" [output_filename.pdf] [jobdir]")
	print("")
	sys.exit(1)

#List all parameters of job
jd = JobData(sys.argv[2])
jd_flat = jd.get_flattened_data()
#for key in jd_flat:
	#print(key, '->', jd_flat[key])		
jd_raw = jd.get_job_raw_data()

output=jd_raw['output']
runtime=jd_raw['jobgeneration']
runtime=runtime['runtime']
print(output)

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
	mmode=int(output[mode+"nmode"])
	ampl=float(output[mode+"ampl"])
	print(i, nmode, mmode, ampl)
	nmodes.append(nmode)
	mmodes.append(mmode)
	ampls.append(ampl)



#Read output_nm files
energy_file    = jd_flat['runtime.p_job_dirpath']+"/output_spec_energy_t00000000000.00000000.txt"
enstrophy_file = jd_flat['runtime.p_job_dirpath']+"/output_spec_enstrophy_t00000000000.00000000.txt"


eps=10e-15
scales = {}
timerescale=1.0
#Remove modes with null values
exclude = ['timestamp']

df_energy=pd.read_csv(energy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
print(df_energy)
df_ens=pd.read_csv(enstrophy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
print(df_ens)

time=df_energy['timestamp']

df_energy_clean = df_energy.loc[:, (df_energy > eps).any(axis=0)]
scales['0']=eps
print(df_energy)
#df_energy.set_index('timestamp',drop=True,inplace=True)
print(df_energy_clean)

exit(1)

##########################################################
# Plotting starts here
##########################################################

print("*"*80)
print("*"*80)
print("*"*80)


fontsize=18
figsize=(10, 10)

fig, axs = plt.subplots(3, figsize=(10,10), sharex=True)
#plt.rc('text', usetex=True)
title="Normal Mode Nonlinear Interaction\n"+params+"\n"+waves
fig.suptitle(title)

for i, ax in enumerate(axs):
	ax.set_xscale("linear", nonposx='clip')
	ax.set_yscale("log", nonposy='clip')
	ylim=[scales[str(i)], 10]
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


c = 0

#ncol=2,handleheight=2.4, labelspacing=0.05
ncol=1
if len(df_geo.columns):
	ncol = 2	
axs[0].set(ylabel='Geostrophic Mode')
df_geo.plot(ax=axs[0])
axs[0].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

ncol=1
if len(df_west.columns):
	ncol = 2	
axs[1].set(ylabel='IGWest Mode')
df_west.plot(ax=axs[1])
axs[1].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

ncol=1
if len(df_east.columns):
	ncol = 2	


df_east.plot(ax=axs[2])
axs[2].set(xlabel="Time (hours)", ylabel='IG East Mode')
axs[2].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

fig.subplots_adjust(right=0.7)

plt.savefig(output_filename, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)

plt.close()


