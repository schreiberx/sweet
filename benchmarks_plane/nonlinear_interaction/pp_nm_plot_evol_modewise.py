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
for key in jd_flat:
	print(key, '->', jd_flat[key])		
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
print("k0:        ", k0)
print("k1:        ", k1)
print("Geo_mode:  ", d0)
print("West_mode: ", dwest)
print("East_mode: ", deast)

params = "h0="+str(runtime['h0'])+", "\
	+"f="+str(runtime['sphere_rotating_coriolis_omega'])+", "\
	+"g="+str(runtime['gravitation'])+", "\
	+"L="+str(runtime['plane_domain_size'])+", "\
	+"M="+str(runtime['space_res_spectral'])




#Read output_nm files
nm_geo_file=jd_flat['runtime.p_job_dirpath']+"/output_nm_geo_evol.txt"
nm_igwest_file=jd_flat['runtime.p_job_dirpath']+"/output_nm_igwest_evol.txt"
nm_igeast_file=jd_flat['runtime.p_job_dirpath']+"/output_nm_igeast_evol.txt"

#Remove modes with null values
exclude = ['n', 'time']

df_geo=pd.read_csv(nm_geo_file, sep='\t', skipinitialspace=True, engine="python")
df_geo=df_geo.loc[:, (df_geo != 0).any(axis=0)]
time=df_geo['time']
df_geo=df_geo.loc[:, df_geo.columns.difference(exclude)]

df_west=pd.read_csv(nm_igwest_file, sep='\t', skipinitialspace=True, engine="python")
df_west=df_west.loc[:, (df_west != 0).any(axis=0)]
df_west=df_west.loc[:, df_west.columns.difference(exclude)]

df_east=pd.read_csv(nm_igeast_file, sep='\t', skipinitialspace=True, engine="python")
df_east=df_east.loc[:, (df_east != 0).any(axis=0)]
df_east=df_east.loc[:, df_east.columns.difference(exclude)]

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

for ax in axs:
	ax.set_xscale("linear", nonposx='clip')
	ax.set_yscale("log", nonposy='clip')


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

axs[0].set(ylabel='Geostrophic Mode')
df_geo.plot(ax=axs[0])
axs[0].legend(loc='center left', bbox_to_anchor= (1.01, 0.5))

axs[1].set(ylabel='IGWest Mode')
df_west.plot(ax=axs[1])
axs[1].legend(loc='center left', bbox_to_anchor= (1.01, 0.5))

axs[2].set(xlabel="Time (seconds)", ylabel='IG East Mode')
df_east.plot(ax=axs[2])
axs[2].legend(loc='center left', bbox_to_anchor= (1.01, 0.5))


plt.savefig(output_filename, transparent=True) #, bbox_inches='tight', pad_inches=0)

plt.close()

plt.show()
exit()




for key, d in data.items():
	x = d['x_values']
	y = d['y_values']
	l = key.replace('_', '\\_')

	print(" + "+l)
	print(x)
	print(y)
	ax.plot(x, y, marker=markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)], label=l)

	c = c + 1


if title != '':
	plt.title(title, fontsize=fontsize)

plt.xlabel("Timestep size $\Delta t$ (sec)", fontsize=fontsize)

#
# Name of data
#
dataname = "TODO"
if 'prog_h' in muletag:
	dataname = "surface height $h$"

#
# Norm
#
if 'linf' in muletag:
	norm = "$L_\infty$"
else:
	norm = "$L_{TODO}$"

plt.ylabel(norm+" error on "+dataname, fontsize=fontsize)



plt.legend(fontsize=15)

plt.savefig(output_filename, transparent=True, bbox_inches='tight', pad_inches=0)

plt.close()

