#
# This program reads in data files from
# of Mark's dns code and makes contour plots
#
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys

from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) < 2:
	print "usage: "+sys.argv[0]+" [infile]"
	sys.exit(1)

in_filename = sys.argv[1]

data=np.fromfile(in_filename,dtype=np.dtype('d'),count=-1,sep="")

#Convert this to meaningful data
TotalBytes = data.nbytes # This is the total size of the data read - in bytes
BytesPerDouble = data.itemsize # This is the total number of bytes per double
TotalItems = TotalBytes/BytesPerDouble
print 'Total Items = ', TotalItems

# The first element in the array is the time at which the snapshot occurs
SnapshotTime = data[0]
print "Time = ", SnapshotTime

# Nx, Ny, Nz
nx = int(data[1])
ny = int(data[2])
nz = int(data[3])
print "nx, ny, nz ", nx, ny, nz

# Declare arrays and store x, y and z coordinates

xcoord = np.zeros(nx)
ycoord = np.zeros(ny)
print "Length of dimensioned xcoord",len(xcoord)
print "Length of dimensioned ycoord",len(ycoord)

if nz != 1:
    zcoord = np.zeros(nz)
    print "Length of dimensioned zcoord",len(zcoord)

xcoord = data[4:3+nx]
ycoord = data[4+nx:3+nx+ny]
if nz != 1:
    zcoord = data[4+nx+ny:3+nx+ny+nz] #For 3D arrays

# Declare arrays and store the data

if nz != 1:
    DataD = np.zeros((nx-1,ny-1,nz-1))
else:
    DataD = np.zeros((nx-1,ny-1))

print "Shape of DataD", DataD.shape
print "DataD elements", DataD[0,0], DataD[1,1]

# Store the data in a way that makes more sense

icounter = 3+nx+ny+nz

if nz !=1:
    for k in np.arange(1,nz):
        for i in np.arange(1,nx):
            for j in np.arange(1,ny):
                icounter += 1
                print k,i,j
                DataD[i-1,j-1,k-1] = data[icounter]
else:
    for i in np.arange(1,nx):
        for j in np.arange(1,ny):
            icounter += 1
            DataD[i-1,j-1]=data[icounter]


# Generate raw binary file
DataD.tofile(in_filename+".raw")

# generate csv file
fd = open(in_filename+".csv", 'w')
for d in DataD:
	fd.write("\t".join([str(i) for i in d]))
	fd.write("\n")
