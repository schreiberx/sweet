#! /usr/bin/env python2

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys
import math


def loadDataFromFile(filename):
	global prefix

	try:
		data = np.loadtxt(filename, skiprows = 0)
	except:
		prefix = filename if len(sys.argv) <= 3 else sys.argv[3]
		print(prefix+": UNABLE TO OPEN "+filename)
		sys.exit(1)

	#labelsx = data[0,1:]
	#labelsy = data[1:,0]
	#data = data[1:,1:]
	return data

ref_data = loadDataFromFile(sys.argv[1])
cmp_data = loadDataFromFile(sys.argv[2])
prefix = sys.argv[3]

norm_l1_value = 0.0
norm_l2_value = 0.0
norm_linf_value = 0.0

#Since data is on centred A-grids, grids with different size do not align perfectly, so are not comparable
# Check grid compatibility
size_j = len(cmp_data)
size_i = len(cmp_data[0])
size_j_ref = len(ref_data)
size_i_ref = len(ref_data[0])
#print "Ref grid:", size_i_ref, "x", size_j_ref, " Computed grid:", size_i, "x", size_j

if size_i == size_i_ref:
        factori=(size_i_ref)/(size_i)
        rem=(size_i_ref)%(size_i)
        #print "Grid Factor i:",factori,"Rem:",rem
        if rem !=0:
                print("Grids not aligned, cannot do comparison")
                sys.exit(1)
else:
        print("Reference has smaller grid than computed data, cannot compate results")
        sys.exit(1)

if size_j <= size_j_ref:
        factorj=(size_j_ref)/(size_j)
        rem=(size_j_ref)%(size_j)
        #print "Grid Factor j:",factorj,"Rem:",rem
        if rem !=0:
                print("Grids not aligned, cannot do comparison")
                sys.exit(1)
elif size_j > size_j_ref:
        print("Reference has smaller grid than computed data, cannot compate results")
        sys.exit(1)

#print cmp_data[1,0], cmp_data[1,1], cmp_data[1,2]
#print ref_data[1,0],  ref_data[1,1],  ref_data[1,2]
#print ref_data[1,0*factori], ref_data[1,1*factori], ref_data[1,2*factori]
for j in range(1, size_j-1):
	for i in range(1, size_i-1):
		value = cmp_data[j,i]-ref_data[j*factorj,i*factori]

		norm_l1_value += abs(value)
		norm_l2_value += value*value
		norm_linf_value = max(norm_linf_value, abs(value))

norm_l1_value = norm_l1_value/(size_i*size_j)
norm_l2_value = math.sqrt(norm_l2_value)/(size_i*size_j)

#
# L1, L2, Linf
#
print(prefix+"\t"+str(norm_l1_value)+"\t"+str(norm_l2_value)+"\t"+str(norm_linf_value))
