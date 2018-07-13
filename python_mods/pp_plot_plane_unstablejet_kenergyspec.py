#! /usr/bin/python3
# Plot unstable jet fields
# 
#--------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import scipy
from scipy.interpolate import RectBivariateSpline

import numpy as np
import sys

#Figure definitions
fontsize=18
figsize=(9, 7)

for filename in sys.argv[1:]:

	#Load data
	print(filename)
	data = np.loadtxt(filename)
	
	if 'spec' not in filename:
		print("Use this only to plot spectrums")
		sys.exit()

	#Get max/min
	cmin = np.amin(data)
	cmax = np.amax(data)
	
	#Set physical grid for axis
	n = data.shape[0]
	m=int(n/2)+1
	#m=100
	x_min = 0
	x_max = int(n/2)
	y_min = 0
	y_max = int(n/2)
	x = np.linspace(x_min, x_max, m)
	y = np.linspace(y_min, y_max, m)
	
	#Labels
	labelsx = np.linspace(x_min, x_max, 10)
	labelsy = np.linspace(y_min, y_max, 10)

	#Set tittle
	title=""	
	if 'ke' in filename:
		title+="Kinetic Energy Spectrum "
		
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
	title += ' days '
	
	if time > 0 :
		title+=" dt="
		pos1 = filename.find('_C')
		pos2 = filename.find('_R')
		title += filename[pos1+2:pos2]
		
	print(title)

	data=data[0:m,0:m]
	
	#Start plotting 2d figure
	plt.figure(1, figsize=figsize)
	
	#2D plot
	datalog=data+1
	datalog=np.log(datalog)
	cmin = np.amin(datalog)
	cmax = np.amax(datalog)
	extent = (labelsx[0], labelsx[-1], labelsy[0], labelsy[-1])
	plt.imshow(datalog, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('jet'))
	plt.clim(cmin, cmax)	
	cbar = plt.colorbar()	
			
	if 'ke' in filename:	
		cbar.set_label('$m^2s^{-2}$', rotation=270, labelpad=+20, size=fontsize)
		cbar.ax.tick_params(labelsize=fontsize) 
	
	
	plt.title(title, fontsize=fontsize)

	#Axis
	ax = plt.gca()
	ax.xaxis.set_label_coords(0.5, -0.075)
	
	plt.xticks(fontsize=fontsize)
	plt.yticks(fontsize=fontsize)
 
	#plt.xticks(labelsx, fontsize=fontsize)
	plt.xlabel("x mode", fontsize=fontsize)

	#plt.yticks(labelsy, fontsize=fontsize)
	plt.ylabel("y mode", fontsize=fontsize)


	#Save file as eps
	outfilename = filename.replace('.csv', '_2D.eps')
	
	outfilename = outfilename.replace('/output', '')
	print(outfilename)
	plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', \
                        pad_inches=0)
	
	plt.close()


	#Start plotting 1d figure
	plt.figure(2, figsize=figsize)
	
	#Polar coordinates for shell spectrum
	
	r=np.arange(0, m, 1) #radius
	t=np.arange(0, np.pi/2, np.pi/2/m) #angle
	
	#Create cubic interpolation of reference file
	interp_spline = RectBivariateSpline(x, y, data)
	
	energy=np.zeros(m)
	i=0
	for ri in r:
		energy[i]=0.0
		j=0
		for tj in t:
			x_tmp=ri*np.cos(tj)
			y_tmp=ri*np.sin(tj)
			energy[i] += interp_spline(x_tmp, y_tmp)/float(m)
			#print(j, x_tmp, y_tmp, energy[i], data[0,0], m)
			j=j+1
		i=i+1
	#Quick check to see if things match
	print("Energy in shells: ", energy[0:10])
	print("Energy in first column of data: ", data[0:10,0])
	print("Shell modes: ", r[0:10])
	
	r_ref53=r[int(m/4)-100:int(m/4)+100]
		
	offsetx=10
	offsety=0.1
	en_ref53=np.array([])
	for tmp in r_ref53:
		ytmp=np.power(offsety*tmp, -float(5.0/3.0))*offsetx
		#print(tmp, ytmp)
		en_ref53=np.append(en_ref53, [ytmp])
		#print(en_ref53)
	
	plt.title(title, fontsize=fontsize)
	
	plt.loglog(r,energy)
	plt.loglog(r_ref53,en_ref53, '--', color='black')
	
	
	#Axis
	ax = plt.gca()
	ax.annotate("-5/3", xy=(float(m/4), en_ref53[100]+0.1), fontsize=fontsize)
	ax.xaxis.set_label_coords(0.5, -0.075)
	ax.set_ylim(0.001, 100)
	#ax.set_yscale('log')

	plt.xticks(fontsize=fontsize)
	plt.yticks(fontsize=fontsize)
 
	#plt.xticks(labelsx, fontsize=fontsize)
	plt.xlabel("horizontal wavenumber", fontsize=fontsize)

	#plt.yticks(labelsy, fontsize=fontsize)
	plt.ylabel("Energy $m^2s^{-2}$", fontsize=fontsize)


	#Save file as eps
	outfilename = filename.replace('.csv', '_1D.eps')
	
	outfilename = outfilename.replace('/output', '')
	print(outfilename)
	plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', \
                        pad_inches=0)
	
	plt.close()

