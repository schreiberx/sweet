#! /usr/bin/python3
# Plot unstable jet fields
# 
#--------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker


import scipy
from scipy.interpolate import RectBivariateSpline

import numpy as np
import sys
from SWEETParameters import *


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
	
	earth = EarthMKSDimensions()
	benchpar = Unstablejet()

	#Domain
	xL_min = benchpar.x_min
	xL_max = benchpar.x_max
	yL_min = benchpar.y_min
	yL_max = benchpar.y_max

	#Set physical grid for axis
	n = data.shape[1]
	#m=int(n/2)+1
	m=int(2*n/3)-1 #anti-aliasing cut
	#m=10
	x_min = 0
	x_max = int(m)
	y_min = 0
	y_max = int(m)
	x = np.linspace(x_min, x_max, m+1)
	y = np.linspace(y_min, y_max, m+1)
	
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


	#Start plotting 2d figure
	plt.figure(1, figsize=figsize)
	
	#2D plot
	datalog=data[0:m,0:m]+1
	datalog=np.log(datalog)
	cmin = np.amin(datalog)
	cmax = np.amax(datalog)
	extent = (labelsx[0], labelsx[-1], labelsy[0], labelsy[-1])
	plt.imshow(datalog, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('jet'))
	plt.clim(cmin, cmax)	
	cbar = plt.colorbar()	
			
	if 'ke' in filename:	
		cbar.set_label('Spectral density ($m^2s^{-2}$)', rotation=270, labelpad=+20, size=fontsize)
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
	
	#Calculate energy per shell
	r=np.arange(0, m+1, 1) #radius
	energy=np.zeros(m+1)
	shell_pattern=np.zeros((m+1, m+1))
	if 0:
		#Polar coordinates for shell spectrum
		t=np.arange(0, np.pi/2, np.pi/2/m) #angle
		
		#Create cubic interpolation of reference file
		interp_spline = RectBivariateSpline(x, y, data)
		
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
	else:
		print("Generating energy in shells (Each x is 1/", m, ")")
		for i in range(0,m):
			for j in range(0,m):
				k=np.sqrt(pow(float(i),2)+pow(float(j),2))
				intk=int(k)
				if intk < m :
					energy[intk]=energy[intk]+data[i,j]
					shell_pattern[i,j]=intk
			print(".", end='', flush=True)
			#print(i, j, k, intk, data[i,j], energy[intk], data.shape, energy.shape)
				
	print(".")
	#Quick check to see if things match
	print("Energy in shells: ", energy[0:10])
	print("Energy in first column of data: ", data[0:10,0])
	print("Shell modes: ", r[0:10])
	print("Pattern:\n", shell_pattern[0:10,0:10])
	
	
	
	#r_ref53=r[-300:-200]
	
	#offsetx=10
	#offsety=0.01
	#en_ref53=np.array([])
	#for tmp in r_ref53:
	#	ytmp=np.power(offsety*tmp, -float(5.0/3.0))*offsetx
		#print(tmp, ytmp)
	#	en_ref53=np.append(en_ref53, [ytmp])
		#print(en_ref53)
	
	r_ref3=r[-int(m/4):-1]
	
	offsetx=m*1000
	offsety=0.005
	en_ref3=np.array([])
	for tmp in r_ref3:
		ytmp=np.power(tmp, -float(5.0/3.0))*offsetx
		#print(tmp, ytmp)
		en_ref3=np.append(en_ref3, [ytmp])
		#print(en_ref53)
	
	plt.title(title, fontsize=fontsize)
	
	#Convert wavenumber to wavelength
	r=xL_max*1000/r
	r_ref3=xL_max*1000/r_ref3
	
	plt.loglog(r, energy)
	#plt.loglog(r_ref53, en_ref53, '-.', color='black')
	plt.loglog(r_ref3, en_ref3, '-.', color='black')
	
	
	#Axis
	ax = plt.gca()
	#ax.annotate("-5/3", xy=(r_ref53[50], en_ref53[50]+0.01), fontsize=fontsize)
	ax.annotate("$k^{-5/3}$", xy=(r_ref3[-int(m/8)], en_ref3[-int(m/8)]+0.001), fontsize=fontsize)

	ax.xaxis.set_label_coords(0.5, -0.075)

	plt.gca().invert_xaxis()

	plt.xticks(fontsize=fontsize)
	plt.yticks(fontsize=fontsize)
 
	#plt.xticks(labelsx, fontsize=fontsize)
	plt.xlabel("Horizontal wavenumber", fontsize=fontsize)
	
	plt.xlabel("Horizontal wavelength (km)", fontsize=fontsize)

	#plt.yticks(labelsy, fontsize=fontsize)
	plt.ylabel("Spectral density $m^3s^{-2}$", fontsize=fontsize)


	#Save file as eps
	outfilename = filename.replace('.csv', '_1D.eps')
	
	outfilename = outfilename.replace('/output', '')
	print(outfilename)
	plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', \
                        pad_inches=0)
	
	plt.close()

