#! /usr/bin/env python3

import sys

import matplotlib
matplotlib.use('Agg')

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker


from mule.postprocessing.JobsData import *

if len(sys.argv) >= 3:
	pickle_tag = sys.argv[1]
	output_filename = sys.argv[2]

else:
	print("")
	print("Usage:")
	print("")
	print("	"+sys.argv[0]+" [pickle tag] [output_filename.pdf] [jobdir1] [jobdir2] ... [jobdirN]")
	print("")
	sys.exit(1)


if len(sys.argv) >= 4:
	# Load Jobs specified via program parameters
	jd = JobsData(job_dirs=sys.argv[3:])
else:
	# Load all Jobs
	jd = JobsData()


data = jd.get_flattened_data()



def label(d):
	val = d['runtime.timestepping_method'].replace('_', '\\_')+', $\Delta t = '+str(d['runtime.timestep_size'])+'$'
	return val

def x_values(d):
	return d[pickle_tag+'.spectrum_wavelength']

def y_values(d):
	return d[pickle_tag+'.spectrum']


##########################################################
# Plotting starts here
##########################################################

print("*"*80)
print("*"*80)
print("*"*80)


fontsize=18
figsize=(9, 7)

plt.figure(1, figsize=figsize)
#plt.rc('text', usetex=True)

ax = plt.gca()

if False:
	# TODO
	# TODO: Plot reference slopes
	# TODO
	en_ref53=np.array([])
	iref=(energy[1]/50.0)/np.power(float(i), -float(5.0/3.0))	
	for tmp in r_ref3:
		ytmp=np.power(tmp, -float(5.0/3.0))*iref
		en_ref53=np.append(en_ref53, [ytmp])

	r_ref3_len=xL_max*1000/r_ref3[1:]
	#plt.loglog(r_ref53, en_ref53, '-.', color='black')
	plt.loglog(r_ref3_len, en_ref3[1:], '-.', color='black')
	plt.loglog(r_ref3_len, en_ref53[1:], '-.', color='black')

	ax.annotate("$k^{-5/3}$", xy=(r_ref3_len[-1]-10, en_ref53[-1]), fontsize=fontsize)
	ax.annotate("$k^{-3}$", xy=(r_ref3_len[-1]-10, en_ref3[-1]), fontsize=fontsize)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_label_coords(0.5, -0.075)
ax.set_facecolor('xkcd:white')

# invert axis for wavelength
ax.invert_xaxis()

ax.xaxis.set_label_coords(0.5, -0.075)

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

#plt.xticks(labelsx, fontsize=fontsize)
plt.xlabel("Horizontal wavelength ($km$)", fontsize=fontsize)

#plt.yticks(labelsy, fontsize=fontsize)
plt.ylabel("H' spectrum", fontsize=fontsize)


markers = ['']
linestyles = ['-', '--', ':', '-.']

c = 0

title = ''
for key, d in data.items():
	print("Processing job in directory: "+d['jobgeneration.job_dirpath'])

	l = label(d)
	x = x_values(d)
	y = y_values(d)

	print(" + "+l)
	plt.loglog(x[1:], y[1:], markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)], label=l)	
	c = c + 1


if title != '':
	plt.title(title, fontsize=fontsize)


plt.legend(fontsize=15)

print("Writing file '"+output_filename+"'")
plt.savefig(output_filename, transparent=True, bbox_inches='tight', pad_inches=0)

plt.close()

