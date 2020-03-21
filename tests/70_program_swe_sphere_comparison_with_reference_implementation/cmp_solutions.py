#! /usr/bin/env python3

import numpy as np
import sys

if len(sys.argv) < 3:
	print("")
	print("Usage:")
	print("	"+argv[0]+" [file1.csv] [file2.csv]")
	print("")

filename_ref = sys.argv[1]
filename_sweet = sys.argv[2]




#filename_ref = 'job_benchref_solution/'+filename
data_ref = np.loadtxt(filename_ref)
#print("Reference data size:")
#print(data_ref.shape)



#filename_sweet = 'job_bench_sweet/'+filename
data_sweet = np.loadtxt(filename_sweet)
# Skip first row and col since they contain the angles
data_sweet = data_sweet[1:,1:]

#print(data_ref[0]-data_sweet[0])
#sys.exit(1)
#print(data_ref-data_sweet)

#print("Sweet data size:")
#print(data_sweet.shape)

diff = data_ref - data_sweet


if 0:
    import matplotlib.pyplot as plt
    plt.imshow(diff)
    plt.colorbar()
    plt.show()
    sys.exit(1)

lmax_error = np.max(np.abs(data_ref-data_sweet))
print(" + Error Lmax: "+str(lmax_error))


if 'prog_phi' in filename_ref:
	if lmax_error > 1e-5:
		raise Exception("ERROR threshold too large!")
else:
	if lmax_error > 1e-8:
		raise Exception("ERROR threshold too large!")

