#! /bin/python

import math
import copy
import numpy as np
import scipy
import sys

print 'Vertical Modes Analysis'

#Max levels
K=4
#K=10
m=0.5
      
lamb=-1.0
 
A=np.zeros((K,K))
for i in range(1,K-1):
	print i
	A[i][i-1]=1
	A[i][i]=-2.*lamb
	A[i][i+1]=1
print "lkjasdflkjasfd"

A=np.delete(A,(0), axis=0)
A=np.delete(A,(0), axis=1)
A=np.delete(A,len(A)-1, axis=0)
A=np.delete(A,len(A[0])-1, axis=1)

print "lkjasdflkjasfd"

e, v=np.linalg.eig(A)

import matplotlib.pyplot as plt
print e
plt.plot(e)
plt.show()
