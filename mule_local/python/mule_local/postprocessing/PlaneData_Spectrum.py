#! /usr/bin/env python3

import sys
import math
import os

import numpy as np

from mule.postprocessing.JobData import *
from mule_local.postprocessing.PlaneDataPhysical import *
from mule.InfoError import *



class PlaneData_Spectrum(InfoError):

    def __init__(
        self,
        filename_phys = None,
        params = []
    ):
        """
        Parameters:

        filename_phys: string
        	Filename (full path) of physical data

        params: list of strings
        	-
        """

        InfoError.__init__(self, "PlaneData_Spectrum")
        self.params = params


        # values of spectrum
        self.spectrum = None

        # mode wavelength
        self.spectrum_wavelength = None

        # mode number
        self.spectrum_mode = None


        if filename_phys != None:
        	self.compute_spectrum(filename_phys, params)


    def compute_spectrum(
        	self,
        	filename_phys,
        	params,
        ):

        #
        # Get meta information about job
        #

        self.job_metadata_available = False

        try:
        	# Step 1) Determine job directory name 
        	jobdir = os.path.dirname(filename_phys)

        	print("Loading data from job directory '"+jobdir+"'")

        	# Step 2) Load job meta information
        	j = JobData(jobdir=jobdir)

        	self.job_metadata_available = True

        except Exception as e:
        	print(str(e))
        	self.job_metadata_available = False


        if self.job_metadata_available:
        	data = j.get_flattened_data()

        	if 'runtime.plane_domain_size' not in data:
        		raise Exception("Physical domain size must be specified via runtime parameters in MULE")

        	domain_size = data['runtime.plane_domain_size']

        	if isinstance(domain_size, list):
        		domain_size = domain_size[0]



        #
        # The following code is converted from a development of Pedro Peixoto to fit into the MULE framework
        #

        # some parameter
        mmin = 0

        # Load file
        udata_ = PlaneDataPhysical(filename_phys)
        udata = udata_.data

        #Calculate spectrum
        #-----------------------------------
        
        print("Physical shape")
        print(udata.shape)

        uspec = np.fft.fft2(udata)/(udata.shape[0]*udata.shape[1])

        print("Spectral shape")
        print(uspec.shape)

        # Calculate amplitude spectrum
        data = np.multiply(uspec,np.conjugate(uspec))
        data = data.real

        n = data.shape[0]

        # Since data u,v data is real, the spectrum has a symmetry and all that matters is the 1st quadrant
        # we multiply by 2 to account for the symmetry
        data = 2*data[0:int(n/2)-1, 0:int(n/2)-1]

        # Adjust data size
        n = data.shape[0]
        
        # m=int(n/2)+1
        m = int(2*n/3)+1 #anti-aliasing cut
        if mmin == 0:
        	mmin = m
        else:
        	if m > mmin:
        		m = mmin

        print("Anti-aliased spectrum region:", m)		

        #Calculate energy per shell
        # TODO: convert to linspace
        r = np.arange(0, m+1, 1) # radius
        energy = np.zeros(m+1)
        shell_pattern = np.zeros((m+1, m+1))

        print("Generating energy in shells (Each . is 1/", m, ")")
        for i in range(0,m):
        	for j in range(0,m):
        		k = np.sqrt(pow(float(i),2)+pow(float(j),2))
        		intk = int(k)
        		if intk < m :
        			energy[intk] = energy[intk]+data[i,j]
        			shell_pattern[i,j] = intk
        	print(".", end='', flush=True)
        	#print(i, j, k, intk, data[i,j], energy[intk], data.shape, energy.shape)

        print(".")	

        #Quick check to see if things match
        #print("Energy in shells: ", energy[0:10])
        #print("Energy in first column of data: ", data[0:10,0])
        #print("Shell modes: ", r[0:10])
        #print("Pattern:\n", shell_pattern[0:10,0:10])

        self.spectrum = energy[:]
        self.spectrum_mode = r[:]

        if self.job_metadata_available:
        	# Convert wavenumber to wavelength
        	self.spectrum_wavelength = np.zeros(m+1)
        	self.spectrum_wavelength[1:] = domain_size*1000/r[1:]

        else:
        	self.spectrum_wavelength = None



    def print(self):

        print("")

        if self.job_metadata_available:
        	for i in range(len(self.spectrum)):
        		print(str(self.spectrum_wavelength[i])+":\t"+str(self.spectrum[i])+"\t"+str(self.spectrum_mode[i]))

        else:
        	for i in range(len(self.spectrum)):
        		print(str(self.spectrum_mode[i])+":\t"+str(self.spectrum[i]))

        print("")



    def write_file(
        	self,
        	picklefile,
        	tagname = None
        ):

        #
        # If picklefile is specified, write norm data to pickle file.
        # This can be later on further postprocessed!
        #
        if picklefile != None:
        	import pickle

        	if tagname != None:
        		tagname += '.'
        	else:
        		tagname = ''

        	pickle_data = {
        		tagname+'spectrum' : self.spectrum,
        		tagname+'spectrum_mode' : self.spectrum_mode,
        		tagname+'spectrum_wavelength' : self.spectrum_wavelength,
        	}

        	print(" + picklefile: "+str(picklefile))

        	with open(picklefile, 'wb') as f:
        		# Pickle the 'data' dictionary using the highest protocol available.
        		pickle.dump(pickle_data, f)

        print("")



if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("")
        print("Usage:")
        print("	"+sys.argv[0]+" [plane_physical_data] [picklefile output file (optional)] [reference tagname (optional)]")
        print("")
        print("	plane_physical_data:")
        print("		filename with physical plane data")
        print("")
        print("	picklefile:")
        print("		filename to write pickle file to")
        print("		spectrum: list")
        print("		spectrum_modes: list")
        print("")
        print(" reference tagname:")
        print("		How to name value in .pickle file")
        print("")
        sys.exit(1)


    ufile = sys.argv[1]

    picklefile = None
    if len(sys.argv) >= 3:
        picklefile = sys.argv[2]

    tagname = None
    if len(sys.argv) >= 4:
        tagname = sys.argv[3]

    s = PlaneData_Spectrum()
    s.compute_spectrum(ufile, sys.argv[3:])

    s.print()
    s.write_file(picklefile, tagname)
