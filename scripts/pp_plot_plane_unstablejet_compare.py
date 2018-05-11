#! /usr/bin/python3
# Plot unstable jet fields
#Campare against reference solutions
#
#--------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys
import subprocess
import os
import math
import scipy
from scipy.interpolate import RectBivariateSpline
import matplotlib.colors as colors


#SWeet stuff
from SWEETParameters import *

#Go to working directory (need because .py is used as link)
retval = os.getcwd()
#print ("Current working directory %s" % retval)
os.chdir( retval )

#Figure definitions
fontsize=16
figsize=(9, 7)

#Arguments check
if len(sys.argv) != 3 :
	print("I need 2 arguments:")
	print("First file is reference file, second file is the analysed file")
	sys.exit(1)


#Get data files
def loadDataFromFile(filename):
	try:
		data = np.loadtxt(filename)
	except:
		print(": UNABLE TO OPEN '"+filename+"'")
		sys.exit(1)

	#labelsx = data[0,0:]
	#labelsy = data[0:,0]
	#data = data[1:,1:]
	return data


filename_ref = sys.argv[1]
data_ref = loadDataFromFile(sys.argv[1])
filename_cmp = sys.argv[2]
data_cmp = loadDataFromFile(sys.argv[2])

#Extract parameters from filename (simplifies analysis)
p_ref = ParameterFilename(filename_ref)
p_cmp = ParameterFilename(filename_cmp)


#Check variable to be plotted
title = p_ref.title
outfile = p_ref.var
varplot = p_ref.var

#Check dimensions
(ny_ref, nx_ref) = data_ref.shape
(ny_cmp, nx_cmp) = data_cmp.shape
#print("Dimensions (ref and cmp)")
#print(ny_ref, nx_ref)
#print(ny_cmp, nx_cmp)

#Set benchmark
earth = EarthMKSDimensions()
benchpar = Unstablejet()

#Domain
x_min = benchpar.x_min
x_max = benchpar.x_max
y_min = benchpar.y_min
y_max = benchpar.y_max

multiplier_j = (ny_ref)/(ny_cmp)
multiplier_i = (nx_ref)/(nx_cmp)

#print ("Multipliers: ", multiplier_i, multiplier_j)

norm_l1_value = 0.0
norm_l2_value = 0.0
norm_linf_value = 0.0

if multiplier_i == 1 and multiplier_j == 1: #Grids are the same
	data = data_cmp - data_ref
elif multiplier_i >= 1 and multiplier_j >= 1 :
	#Comparison via interpolation
	#print("Interpolation")
	# A-grid REFERENCE (file1) - sweet outputs only A grids physical space
	dx_ref=(x_max-x_min)/(nx_ref)
	dy_ref=(y_max-y_min)/(ny_ref)
	x_ref = np.arange(x_min, x_max, dx_ref)
	y_ref = np.arange(y_min, y_max, dy_ref)
	x_ref += dx_ref/2
	y_ref += dy_ref/2
	X_ref, Y_ref = np.meshgrid(x_ref, y_ref)

	#Create cubic interpolation of reference file
	interp_spline = RectBivariateSpline(y_ref, x_ref, data_ref)

	#A-grid cmp file (file2)
	dx_cmp=(x_max-x_min)/nx_cmp
	dy_cmp=(y_max-y_min)/ny_cmp
	x_cmp = np.arange(x_min, x_max, dx_cmp)
	y_cmp = np.arange(y_min, y_max, dy_cmp)
	x_cmp += dx_cmp/2
	y_cmp += dy_cmp/2
	X_cmp, Y_cmp = np.meshgrid(x_cmp, y_cmp)

	#Get reduced reference resolution
	data_ref_low = interp_spline(y_cmp, x_cmp)
	data = data_cmp - data_ref_low

else :
	print ("Please provide reference solution (file1) with dimension larger or equal file2")
	sys.exit(1)


#norm_l1_value = data
norm_l1_value = np.sum(np.absolute(data))
norm_l1_value = norm_l1_value/(nx_cmp*ny_cmp)
norm_l2_value = np.sum(np.square(data))
norm_l2_value = math.sqrt(norm_l2_value/(nx_cmp*ny_cmp))
norm_linf_value = np.max(np.absolute(data))


#Set physical grid for axis
n = data.shape[0]
x = np.linspace(x_min, x_max, n)
y = np.linspace(y_min, y_max, n)

#Labels
labelsx = np.linspace(x_min, x_max, 7)
labelsy = np.linspace(y_min, y_max, 7)


#Start plotting
plt.figure(figsize=figsize)

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


#Get max/min
cmin = np.amin(data)
cmax = np.amax(data)
mid_val = 0

dx_cmp = (x_max-x_min)/nx_cmp
dy_cmp = (y_max-y_min)/ny_cmp
extent = (labelsx[0]+dx_cmp/2, labelsx[-1]-dx_cmp/2, labelsy[0]+dy_cmp/2, labelsy[-1]-dy_cmp/2)

#Color plot
#plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto')
#plt.imshow(data, interpolation='nearest', origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))
#plt.imshow(data, interpolation='nearest', origin='lower')
plt.imshow(data, interpolation='nearest', extent=extent, origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'), norm=MidpointNormalize(midpoint=mid_val,vmin=cmin, vmax=cmax))
#plt.imshow(ras, cmap=cmap, clim=(elev_min, elev_max), norm=MidpointNormalize(midpoint=mid_val,vmin=elev_min, vmax=elev_max))
#plt.pcolor(x_cmp, y_cmp, data)

#Colorbar - centred
plt.clim(cmin, cmax)
cref=max(abs(cmin),abs(cmax))
plt.clim(-cref, +cref)
if varplot == 'vort':
	cbar = plt.colorbar(format='%.0e')
	cbar.set_label('1/s', rotation=270, labelpad=+20, size=fontsize)
elif varplot == "depth":
	cbar = plt.colorbar()
	cbar.set_label('meters', rotation=270,labelpad=+5, size=fontsize)
else:
	cbar = plt.colorbar()

cbar.ax.tick_params(labelsize=fontsize)


#Contour lines (black)
if varplot == 'vort':
	s = 2e-5
	eta_contour_levels = np.append(np.arange(-1e-4, 0, s), np.arange(s, 1e-4, s))

	#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
	#plt.contour(x,y,data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=eta_contour_levels, linewidths=0.5)
	#plt.contourf(x, y, data, vmin=cmin, vmax=cmax, levels=eta_contour_levels)
elif varplot == "depth":
	#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
	hs = 20
	#h_contour_levels = np.append(np.arange(900, 1000-hs, hs), np.arange(1000+hs, 1100, hs))
	h_contour_levels = np.append(np.arange(-20, 0, hs), np.arange(0, 20, hs))
	#plt.contour(x,y, data, colors="black", origin='lower', vmin=cmin, vmax=cmax, levels=h_contour_levels, linewidths=0.5)
else:
	if cmin != cmax:
		pass
		#plt.contour(data, colors="black", origin='lower', extent=extent, vmin=cmin, vmax=cmax, linewidths=0.5)


#Labels
if p_ref.method == p_cmp.method:
	title+=p_ref.method_paper
	outfile += "_"+p_ref.method
else:
	title += p_ref.method_paper
	title += " vs "
	title += p_cmp.method_paper
	outfile += "_"+p_ref.method
	outfile += "_vs_"
	outfile += p_cmp.method

#Time
title += '  t='
title += str(p_ref.time)
outfile += "_t"
outfile += str(p_ref.time)
if p_ref.time != p_cmp.time:
	title += " vs t="
	title += str(p_cmp.time)
	outfile += "_vs_t"
	outfile += str(p_cmp.time)
title += ' days '

#Time step
title += " dt="
outfile += "_dt"
outfile += str(p_ref.timestep)
if p_ref.method_paper != "REF":
	title += str(p_ref.timestep)
	
if p_ref.timestep != p_cmp.timestep:
	if p_ref.method_paper != "REF":
		title += " vs "
	title += p_cmp.timestep
	outfile += "_vs_dt"
	outfile += str(p_cmp.timestep)
title += ' sec '


print(title)
plt.title(title, fontsize=fontsize)
#plt.title(filename, fontsize=fontsize)

#Axis
ax = plt.gca()
ax.xaxis.set_label_coords(0.5, -0.075)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

#plt.xticks(labelsx, fontsize=fontsize)
plt.xlabel("x (1000 km)", fontsize=fontsize)

#plt.yticks(labelsy, fontsize=fontsize)
plt.ylabel("y (1000 km)", fontsize=fontsize)

plt.show()

#plt.close()

#Save file as eps
#outfilename = filename.replace('.csv', 'compare.eps')
#print(outfilename)
outfilename = outfile+".eps"
print(outfilename)
#outfilename="tmp.eps"
#print(outfilename)
plt.savefig(outfilename, dpi=300, transparent=True, bbox_inches='tight', \
                        pad_inches=0)


#
# L1, L2, Linf
outfile_errors = outfilename.replace('.eps', 'Errors12max.txt')
print(outfile_errors+"\t"+str(norm_l1_value)+"\t"+str(norm_l2_value)+"\t"+str(norm_linf_value))
print(str(norm_l1_value)+"\t"+str(norm_l2_value)+"\t"+str(norm_linf_value), file=open(outfile_errors, 'w'))
