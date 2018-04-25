#! /usr/bin/env python3

import numpy as np
import sys
import math
import scipy
from scipy.interpolate import RectBivariateSpline

sys.stdout.write('')

#SWeet stuff
from SWEETParameters import *

#Go to working directory (need because .py is used as link)
retval = os.getcwd()
#print ("Current working directory %s" % retval)
os.chdir( retval )


#Arguments check
if len(sys.argv) < 3 :
	print("I need 2 arguments:")
	print("1) Reference file;")
	print("2) Analysed file")
	print("3) Prefix (name of directory of analysed file) - optional")
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

if len(sys.argv) == 4:
	prefix = sys.argv[3]
else:
	prefix = p_ref.details+p_cmp.details

#Check variable to be analysed
if p_ref.var != p_cmp.var:
	print("Warning: you are comparing apples and oranges!!!!")
	print(p_ref.var+" vs "+p_cmp.var)

#Check dimensions
(ny_ref, nx_ref) = data_ref.shape
(ny_cmp, nx_cmp) = data_cmp.shape
#print("Ref: ", ny_ref, "x", nx_ref, " ; Data", ny_cmp, "x", nx_cmp)

multiplier_j = (ny_ref)/(ny_cmp)
multiplier_i = (nx_ref)/(nx_cmp)

norm_l1_value = 0.0
norm_l2_value = 0.0
norm_linf_value = 0.0

if multiplier_i == 1 and multiplier_j == 1: #Grids are the same
	data = data_cmp - data_ref
elif multiplier_i > 1 and multiplier_j > 1 :
	#Comparison via interpolation
		#print("Interpolation")
		# A-grid REFERENCE (file1) - sweet outputs only A grids physical space
		dx_ref=1.0/(nx_ref)
		dy_ref=1.0/(ny_ref)
		x_ref = np.arange(0, 1, dx_ref)
		y_ref = np.arange(0, 1, dy_ref)
		x_ref += dx_ref/2
		y_ref += dy_ref/2
		X_ref, Y_ref = np.meshgrid(x_ref, y_ref)

		#Creat cubic interpolation of reference file
		interp_spline = RectBivariateSpline(y_ref, x_ref, data_ref)

		#A-grid cmp file (file2)
		dx_cmp=1.0/nx_cmp
		dy_cmp=1.0/ny_cmp
		x_cmp = np.arange(0, 1, dx_cmp)
		y_cmp = np.arange(0, 1, dy_cmp)
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

#
# L1, L2, Linf
#
print(prefix+"\t"+str(norm_l1_value)+"\t"+str(norm_l2_value)+"\t"+str(norm_linf_value))
