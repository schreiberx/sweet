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

nm_name=output['simulation_benchmark_normal_modes.case']
nwaves=int(output['simulation_benchmark_normal_modes.nwaves'])

k0=[]
k1=[]
d0=[]
dwest=[]
deast=[]
waves = "Initial waves (k0,k1):(geo,west,east) "
for i in range(nwaves):
	wave="simulation_benchmark_normal_modes.w"+str(i)+"."
	k0.append(int(output[wave+"k0"]))
	k1.append(int(output[wave+"k1"]))
	d0.append(int(output[wave+"d0"]))
	dwest.append(int(output[wave+"dwest"]))
	deast.append(int(output[wave+"deast"]))
	waves = waves+"\n ("+str(k0[i])+","+str(k1[i])+"):("+str(d0[i])+","+str(dwest[i])+","+str(deast[i])+") "

print("Initial conditions")
print(waves)
#print("k1:        ", k1)
#print("Geo_mode:  ", d0)
#print("West_mode: ", dwest)
#print("East_mode: ", deast)

params = "h0="+str(runtime['h0'])+", "\
	+"f="+str(runtime['sphere_rotating_coriolis_omega'])+", "\
	+"g="+str(runtime['gravitation'])+", "\
	+"L="+str(runtime['plane_domain_size'])+", "\
	+"M="+str(runtime['space_res_spectral'])
	#+"N="+str(runtime['space_res_physical'])





#Read output_nm files
nm_geo_file    = jd_flat['runtime.p_job_dirpath']+"/output_nm_geo_evol.txt"
nm_igwest_file = jd_flat['runtime.p_job_dirpath']+"/output_nm_igwest_evol.txt"
nm_igeast_file = jd_flat['runtime.p_job_dirpath']+"/output_nm_igeast_evol.txt"

eps=10e-4
eps2=10e-4
scales = {}
timerescale=60*60
#Remove modes with null values
exclude = ['n']

df_geo=pd.read_csv(nm_geo_file, sep='\t', skipinitialspace=True, engine="python")
time=df_geo['time']
df_geo=df_geo.loc[:, df_geo.columns.difference(exclude)]
df_geo_tmp=df_geo.loc[:, (df_geo > eps).any(axis=0)]
scales['0']=eps
if len(df_geo_tmp.columns) < 2:
	print("Tolerance too large for geo, empty plot! Reducing tolerance!")
	scales['0']=eps*eps2
	df_geo=df_geo.loc[:, (df_geo > eps*eps2).any(axis=0)]
else:
	df_geo = df_geo_tmp
df_geo.time=df_geo.time/timerescale
df_geo.set_index('time',drop=True,inplace=True)

df_west=pd.read_csv(nm_igwest_file, sep='\t', skipinitialspace=True, engine="python")
df_west=df_west.loc[:, df_west.columns.difference(exclude)]
df_west_tmp=df_west.loc[:, (df_west > eps).any(axis=0)]
scales['1']=eps
if len(df_west_tmp.columns) < 2:
	print("Tolerance too large for west, empty plot! Reducing tolerance!")
	scales['1']=eps*eps2
	df_west=df_west.loc[:, (df_west > eps*eps2).any(axis=0)]
else:
	df_west=df_west_tmp
df_west.time=df_west.time/timerescale
df_west.set_index('time',drop=True,inplace=True)

df_east=pd.read_csv(nm_igeast_file, sep='\t', skipinitialspace=True, engine="python")
df_east=df_east.loc[:, df_east.columns.difference(exclude)]
df_east_tmp=df_east.loc[:, (df_east > eps).any(axis=0)]
scales['2']=eps
if len(df_east_tmp.columns) < 2:
	print("Tolerance too large for east, empty plot! Reducing tolerance!")
	scales['2']=eps*eps2
	df_east=df_east.loc[:, (df_east > eps*eps2).any(axis=0)]
else:
	df_east=df_east_tmp
df_east.time=df_east.time/timerescale
df_east.set_index('time',drop=True,inplace=True)

print(df_geo)
print(df_west)
print(df_east)


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


