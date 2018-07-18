#! /usr/bin/python3
# Plot unstable jet fields
# 
#--------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

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

plt.figure(figsize=figsize)


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

for filename in sys.argv[1:]:
	
	print("")
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
	n = data.shape[1]
	#m=int(n/2)+1
	m=int(2*n/3)-1 #anti-aliasing cut
	if mmin == 0:
		mmin = m
	else:
		if m > mmin:
			m=mmin
			
	#m=10

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
			
	pos1 = filename.find('_u')
	pos2 = filename.find('_U')
	dif = filename[pos1+2:pos2]
	print(dif)
	dif = float(dif)
	print("Difusion:", dif)
	
	if dif>0:
		dif=str(dif/1000000)
		method1 += "v"+dif
		
	pos1 = filename.find('output')
	name = filename[pos1:]
	pos2 = name.find('_t')
	pos3 = filename.find('.csv')
	time = filename[pos1+pos2+2:pos3]
	time = float(time)
	time = time / 86400
	if timeold != str(time):
		title += 't='
		timeold = str(time)
		title += str(time)
		title += ' days '
		outfilename += str(time)
		outfilename += 'days'
	
	if time > 0 :
		pos1 = filename.find('_C')
		pos2 = filename.find('_R')
		method1 += "dt"+filename[pos1+2:pos2]
		
	
	print(method1)
	outfilename += str(method1)
	
	#Calculate energy per shell
	r=np.arange(0, m+1, 1) #radius
	energy=np.zeros(m+1)
	shell_pattern=np.zeros((m+1, m+1))
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
	#print("Pattern:\n", shell_pattern[0:10,0:10])
	
	#Convert wavenumber to wavelength
	r=xL_max*1000/r
	
	plt.loglog(r, energy, markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)], label=method1)
	#plt.loglog(r_ref53, en_ref53, '-.', color='black')
	
	c = c + 1
	
#Reference line
r_ref3=r[-int(m/4):-1]

offsetx=m*100000

en_ref3=np.array([])
for tmp in r_ref3:
	ytmp=np.power(tmp, -float(5.0/3.0))*offsetx
	#print(tmp, ytmp)
	en_ref3=np.append(en_ref3, [ytmp])
	#print(en_ref53)	

r_ref3=xL_max*1000/r_ref3
plt.loglog(r_ref3, en_ref3, '-', color='black')

print(title)
plt.title(title, fontsize=fontsize)

#Axis
ax = plt.gca()
ax.set_facecolor('xkcd:white')
#ax.annotate("-5/3", xy=(r_ref53[50], en_ref53[50]+0.01), fontsize=fontsize)
ax.annotate("$k^{-5/3}$", xy=(r_ref3[-int(m/8)], en_ref3[-int(m/8)]+0.001), fontsize=fontsize)

ax.xaxis.set_label_coords(0.5, -0.075)

plt.gca().invert_xaxis()

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
 
#plt.xticks(labelsx, fontsize=fontsize)
#plt.xlabel("Horizontal wavenumber", fontsize=fontsize)
	
plt.xlabel("Horizontal wavelength (km)", fontsize=fontsize)

#plt.yticks(labelsy, fontsize=fontsize)
plt.ylabel("Spectral density ($m^3s^{-2}$)", fontsize=fontsize)

plt.legend(fontsize=15)


#Save file as eps
outfilename += '.eps'
	
print(outfilename)
plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', \
                       pad_inches=0)
	
plt.close()

