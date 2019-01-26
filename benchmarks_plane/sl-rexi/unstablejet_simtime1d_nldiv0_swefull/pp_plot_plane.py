#! /usr/bin/env python3
#
# Plot unstable jet fields
# 
#--------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from mule.postprocessing.JobData import *

import numpy as np
import sys
import os

#Figure definitions
fontsize=18
figsize=(9, 7)

if len(sys.argv) <= 1:
	print("")
	print("Usage:")
	print("")
	print(" $ "+sys.argv[0]+" [physdata1] [physdata2] [...]")
	print("")
	sys.exit(1)

for filename in sys.argv[1:]:

	jd = JobData(os.path.split(filename)[0])
	jobdata = jd.get_flattened_data()

	#Load data
	print(filename)
	data = np.loadtxt(filename)

	if 'prog_h' in filename:
		depth = jobdata['runtime.h0']

		print("Summing full depth")
		print(depth)

		data = data + depth

		#Define in km
		data = data / 1000

	#Get max/min
	cmin = np.amin(data)
	cmax = np.amax(data)
	
	#Set physical grid for axis
	x_min = 0
	x_max = 4.00316e+07/1000/1000
	y_min = 0
	y_max = 4.00316e+07/1000/1000
	n = data.shape[0]
	x = np.linspace(x_min, x_max, n)
	y = np.linspace(y_min, y_max, n)
	
	#Labels
	labelsx = np.linspace(x_min, x_max, 7)
	labelsy = np.linspace(y_min, y_max, 7)

	#Start plotting
	plt.figure(figsize=figsize)
	
	#Contour levels for fields
	extent = (labelsx[0], labelsx[-1], labelsy[0], labelsy[-1])
	
	#Color plot
	#plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto')
	#plt.imshow(data, interpolation='nearest', origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))
	#plt.imshow(data, interpolation='nearest', origin='lower')

	#Colorbar
	if 'diag_vort' in filename:
		plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))

		#Fix max and min for vorticity
		plt.clim(cmin, cmax)
		cmin = -5e-5
		cmax = 5e-5
		s = 2e-5
		eta_contour_levels = np.append(np.arange(-1e-4, 0, s), np.arange(s, 1e-4, s))
		plt.clim(cmin, cmax)
		cref=max(abs(cmin),abs(cmax))
		plt.clim(-cref, +cref)
		cbar = plt.colorbar(format='%.0e')
	elif 'prog_h' in filename:
		plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('jet'))
		plt.clim(cmin, cmax)
		hs = 1000000
		h_contour_levels = np.append(np.arange(900, 1000-hs, hs), np.arange(1000+hs, 1100, hs))
		cbar = plt.colorbar()
	elif 'spec' in filename:
		data=data+1
		data=np.log(data)
		plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('jet'))
		cmin = np.amin(data)
		cmax = np.amax(data)
		plt.clim(cmin, cmax)	
		cbar = plt.colorbar()	
	else:
		plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))
		plt.clim(cmin, cmax)
		cref=max(abs(cmin),abs(cmax))
		plt.clim(-cref, +cref)
		cbar = plt.colorbar()
			
	#Contour lines (black)
	if 'diag_vort' in filename:
		pass
#		plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
		#plt.contour(x,y,data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
		#plt.contourf(x, y, data, vmin=cmin, vmax=cmax, levels=eta_contour_levels)
	elif 'prog_h' in filename:
		pass
		#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
		#plt.contour(x,y, data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
	else:
		if cmin != cmax:
			pass
			#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, linewidths=0.5)

	#Set tittle
	title=""
	if 'diag_vort' in filename:
		title += "Vorticity (1/t)"
		cbar.set_label('1/s', rotation=270, labelpad=+20, size=fontsize)
	if 'prog_h' in filename:
		title += "Depth (km) "
		cbar.set_label('km', rotation=270, size=fontsize)
	if 'prog_u' in filename:
		title += "U-Velocity (m/s) "
		cbar.set_label('m/s', rotation=270, size=fontsize)
	if 'prog_v' in filename:
		title += "V-Velocity (m/s) "
		cbar.set_label('m/s', rotation=270, size=fontsize)
		
	cbar.ax.tick_params(labelsize=fontsize) 
	
	#Method
	print("Methods")
	pos1 = filename.find('_tsm_')
	pos2 = filename.find('_tso')
	method1 = filename[pos1+5:pos2]
	print(method1)

	if method1 == "l_cn_na_sl_nd_settls":
		method1 = "SL-SI-SETTLS"
	elif method1 == "l_rexi_na_sl_nd_settls":
		method1 = "SL-EXP-SETTLS"
	elif method1 == "l_rexi_na_sl_nd_etdrk":
		method1 = "SL-ETD2RK"
	elif method1 == "l_rexi_n_etdrk":
		method1 = "ETD2RK"
	elif method1 == "ln_erk":
		if 'ref' in filename:
			method1 = "REF"
		else:
			method1 = "RK-FDC"
			
	title+= method1
	title+= " "

	title += 't='
	pos1 = filename.find('output')
	name = filename[pos1:]
	pos2 = name.find('_t')
	pos3 = filename.find('.csv')
	time = filename[pos1+pos2+2:pos3]
	time = float(time)
	time = time / 86400
	
	title += str(time)
	title += 'days '
	
	if time > 0 :
	    title += " dt="+str(jobdata['runtime.timestep_size'])

	print(title)
	plt.title(title, fontsize=fontsize)

	#Axis
	ax = plt.gca()
	ax.xaxis.set_label_coords(0.5, -0.075)
	
	plt.xticks(fontsize=fontsize)
	plt.yticks(fontsize=fontsize)
 
	#plt.xticks(labelsx, fontsize=fontsize)
	plt.xlabel("x (1000 km)", fontsize=fontsize)

	#plt.yticks(labelsy, fontsize=fontsize)
	plt.ylabel("y (1000 km)", fontsize=fontsize)
	
	# Save file as eps
	outfilename = filename.replace('.csv', '.pdf')

	outfilename = outfilename.replace('output', 'plot')
	print(outfilename)
	plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
	
	plt.close()
