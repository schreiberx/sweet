#! /usr/bin/python3
# Plot unstable jet fields
# 
#--------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys

#Figure definitions
fontsize=12
figsize=(9, 7)


i=0
for filename in sys.argv[1:]:
	i = i + 1
	print(filename)
	if i == 1:
		filename1 = filename
		data1 = np.loadtxt(filename) 
	else:
		filename2 = filename
		data2 = np.loadtxt(filename) 
	
data = data1 - data2


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
s = 2e-5
eta_contour_levels = np.append(np.arange(-1e-4, 0, s), np.arange(s, 1e-4, s))
hs = 5
h_contour_levels = np.append(np.arange(900, 1000-hs, hs), np.arange(1000+hs, 1100, hs))

extent = (labelsx[0], labelsx[-1], labelsy[0], labelsy[-1])

#Color plot
#plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto')
#plt.imshow(data, interpolation='nearest', origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))
#plt.imshow(data, interpolation='nearest', origin='lower')
plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))

#Colorbar
plt.clim(cmin, cmax)
if 'diag_vort' in filename:
	cbar = plt.colorbar(format='%.0e')
else:
	cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=fontsize) 

#Contour lines (black)
if 'diag_vort' in filename:
#		plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
	plt.contour(x,y,data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
	#plt.contourf(x, y, data, vmin=cmin, vmax=cmax, levels=eta_contour_levels)
elif 'prog_h' in filename:
	#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
	plt.contour(x,y, data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
else:
	if cmin != cmax:
		pass
		#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, linewidths=0.5)

#Set tittle
title=""
if 'diag_vort' in filename:
	title+="Vorticity Deviation"
if 'prog_h' in filename:
	title+="Depth Deviation (m) "
if 'l_cn_na_sl_nd_settls' in filename1 and 'l_cn_na_sl_nd_settls' in filename2:
	title+=" SL-SI-SETTLS "
elif 'l_rexi_na_sl_nd_settls' in filename1 and 'l_rexi_na_sl_nd_settls' in filename2:
	title+="SL-EXP-SETTLS "
else:
	title+="SL-EXP-SETTLS vs SL-SI-SETTLS"
	
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

title+=" dt="
pos1 = filename.find('_C')
pos2 = filename.find('_R')
title += filename[pos1+2:pos2]

print(title)
plt.title(title, fontsize=fontsize)
#plt.title(filename, fontsize=fontsize)

#Axis
ax = plt.gca()
ax.xaxis.set_label_coords(0.5, -0.075)

#plt.xticks(labelsx, fontsize=fontsize)
plt.xlabel("x (1000 km)", fontsize=fontsize)

#plt.yticks(labelsy, fontsize=fontsize)
plt.ylabel("y (1000 km)", fontsize=fontsize)

plt.show()

#Save faile as eps
outfilename = filename.replace('.csv', 'compare.eps')
print(outfilename)
plt.savefig(outfilename, dpi=300)
#plt.show()
plt.close()
	
	
