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

earth = EarthMKSDimensions()
benchpar = Unstablejet()

#Domain
xL_min = benchpar.x_min
xL_max = benchpar.x_max
yL_min = benchpar.y_min
yL_max = benchpar.y_max

timeold=""

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

#markers = []
#for m in Line2D.markers:
#    try:
#        if len(m) == 1 and m != ' ' and m != '':
#            markers.append(m)
#    except TypeError:
#        pass

linestyles = ['-', '--', ':', '-.']
markers = ['.', ',', '*', '+', 'x']
markers = ['']
#['.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_', 'P', 'X']

c = 0
outfilename = "kespectrum_"
mmin=0


if len(sys.argv) <= 1:
	print("Arguments must be files with the zonal velocities")
	sys.exit()
    
    
#Labels
title="Kinetic Energy Spectrum "
	
plt.figure(1, figsize=figsize)

for filename in sys.argv[1:]:
	#Load data
	#-----------------------------------
	
	ufile=filename
	if 'prog_u' not in ufile:
		print("Arguments must be the zonal velocities")
		sys.exit()
	
	vfile=filename.replace('prog_u', 'prog_v')
	
	udata = np.loadtxt(ufile)
	vdata = np.loadtxt(vfile)
	
	#Calculate spectrum
	#-----------------------------------
	
	print("Dimensions (u,v)")
	print(udata.shape, vdata.shape)

	uspec=np.fft.fft2(udata)/udata.shape[0]/udata.shape[1]
	vspec=np.fft.fft2(vdata)/vdata.shape[0]/vdata.shape[1]

	print("Spectral Dimensions (u,v)")
	print(uspec.shape, vspec.shape)

	#Calculate full KE spectrum
	#see https://journals.ametsoc.org/doi/10.1175/1520-0469%282001%29058<0329%3ATHKESA>2.0.CO%3B2
	data=0.5*(np.multiply(uspec,np.conjugate(uspec)))+0.5*(np.multiply(vspec,np.conjugate(vspec)))
	data=data.real

	n=data.shape[0]
	#print(data.shape)

	#Since data u,v data is real, the spectrum has a symmetry and all that matters is the 1st quadrant
	#we multiply by 2 to account for the symmetry
	data=2*data[0:int(n/2)-1, 0:int(n/2)-1]

	#Adjust data size
	n=data.shape[0]
	
	#m=int(n/2)+1
	m=int(2*n/3)+1 #anti-aliasing cut
	if mmin == 0:
		mmin = m
	else:
		if m > mmin:
			m=mmin
	print("Anti-aliased spectrum region:", m)		
	#Naming
	#-----------------------------------------
	
	#Set tittle
	title="Kinetic Energy Spectrum "
			
	#Method
	print("")
	print("Analysisng the method/data:")
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
			
	pos1 = filename.find('_u')
	pos2 = filename.find('_U')
	dif = filename[pos1+2:pos2]
	print(dif)
	dif = float(dif)
	print("Difusion:", dif)
	
	if dif>0:
		dif=str(dif/1000000)
		method1 += "_v"+dif
		
	pos1 = filename.find('output')
	name = filename[pos1:]
	pos2 = name.find('_t')
	pos3 = filename.find('.csv')
	time = filename[pos1+pos2+2:pos3]
	time = float(time)
	time = time / 86400
	if timeold != str(time):
		title += '_t='
		timeold = str(time)
		title += str(time)
		title += ' days '
		outfilename += str(time)
		outfilename += 'days'
	
	if time > 0 :
		pos1 = filename.find('_C')
		pos2 = filename.find('_R')
		method1 += "_dt"+filename[pos1+2:pos2]
		
	if len(sys.argv) == 2:
		title+=method1
		
	print(method1)
	outfilename += str(method1)
	
	if len(sys.argv) == 2:
		
		#2D spectrum plot
		#------------------------------------
		print("")
		print("------------------------------------------------------------------------")
		print("        2D spectrum ")
		#Start plotting 2d figure
		plt.figure(2, figsize=figsize)
			
		x_min = 0
		x_max = int(m)
		y_min = 0
		y_max = int(m)
		x = np.linspace(x_min, x_max, m+1)
		y = np.linspace(y_min, y_max, m+1)

		#Labels
		labelsx = np.linspace(x_min, x_max, 10)
		labelsy = np.linspace(y_min, y_max, 10)

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
			cbar.set_label('Kinetic Energy ($m^2s^{-2}$)', rotation=270, labelpad=+20, size=fontsize)
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
		outfilename2d = filename.replace('.csv', '_2D.eps')

		outfilename2d = outfilename2d.replace('/output', '')
		print(outfilename2d)
		plt.savefig(outfilename2d, dpi=300, transparent=True, bbox_inches='tight', \
							pad_inches=0)

		plt.close()

	print("")
	print("------------------------------------------------------------------------")
	print("        1D spectrum ")
	#Start plotting 1d figure
	plt.figure(1, figsize=figsize)

	#Calculate energy per shell
	r=np.arange(0, m+1, 1) #radius
	energy=np.zeros(m+1)
	shell_pattern=np.zeros((m+1, m+1))
	#print("Generating energy in shells (Each x is 1/", m, ")")
	for i in range(0,m):
		for j in range(0,m):
			k=np.sqrt(pow(float(i),2)+pow(float(j),2))
			intk=int(k)
			if intk < m :
				energy[intk]=energy[intk]+data[i,j]
				shell_pattern[i,j]=intk
		#print(".", end='', flush=True)
		#print(i, j, k, intk, data[i,j], energy[intk], data.shape, energy.shape)
				
	#Quick check to see if things match
	#print("Energy in shells: ", energy[0:10])
	#print("Energy in first column of data: ", data[0:10,0])
	#print("Shell modes: ", r[0:10])
	#print("Pattern:\n", shell_pattern[0:10,0:10])

	#Convert wavenumber to wavelength
	rlen=xL_max*1000/r[1:]

	plt.loglog(rlen, energy[1:], markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)], label=method1)
	c = c + 1
	#plt.loglog(rlen, energy[1:])



#Define reference lines -3 and -5/3
r_ref3=r[-int(m/2):-1]
offsetx=m*1000
offsety=0.005
en_ref3=np.array([])
i=int(r_ref3[0])
iref=(energy[1]/10.0)/np.power(float(i), -3)
for tmp in r_ref3:
	ytmp=np.power(tmp, -float(3.0))*iref
	en_ref3=np.append(en_ref3, [ytmp])

en_ref53=np.array([])
iref=(energy[1]/10.0)/np.power(float(i), -float(5.0/3.0))	
for tmp in r_ref3:
	ytmp=np.power(tmp, -float(5.0/3.0))*iref
	en_ref53=np.append(en_ref53, [ytmp])

r_ref3_len=xL_max*1000/r_ref3[1:]	
#plt.loglog(r_ref53, en_ref53, '-.', color='black')
plt.loglog(r_ref3_len, en_ref3[1:], '-.', color='black')
plt.loglog(r_ref3_len, en_ref53[1:], '-.', color='black')

ax = plt.gca()
ax.annotate("$k^{-5/3}$", xy=(r_ref3_len[-1]-10, en_ref53[-1]), fontsize=fontsize)
ax.annotate("$k^{-3}$", xy=(r_ref3_len[-1]-10, en_ref3[-1]), fontsize=fontsize)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_label_coords(0.5, -0.075)
ax.set_facecolor('xkcd:white')

#invert axis for wavelength
plt.gca().invert_xaxis()

#Sort out labels
#-----------------------------
plt.title(title, fontsize=fontsize, y=1.02)

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

#plt.xlabel("Horizontal wavenumber", fontsize=fontsize)
plt.xlabel("Horizontal wavelength (km)", fontsize=fontsize)

#plt.yticks(labelsy, fontsize=fontsize)
plt.ylabel("Kinetic Energy Spectrum $(m^2s^{-2})$", fontsize=fontsize)
plt.legend(fontsize=15)


#Save file as eps
outfilename += '.eps'
print(outfilename)
plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', \
					pad_inches=0)

plt.close()

