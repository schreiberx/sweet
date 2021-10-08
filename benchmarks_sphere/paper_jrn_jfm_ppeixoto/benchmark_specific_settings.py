# ---------------------------------------------
# Class to setup spherical modes initialization 
# author: Pedro Peixoto <ppeixoto@usp.br>
# Oct 2021
# ----------------------------------------
import numpy as np
import pickle
from numpy.lib.function_base import append
import pandas as pd
import re
import os
import os.path

import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import matplotlib.ticker as mtick

from mule.postprocessing.JobData import *

#-----------------------------------------------------------------
# Test Cases for different settings of modes initialization
#   > Used for pre-processing, to creat jobs for sweet
#-----------------------------------------------------------------

class modes_TC1: #Init with energy in full shells from n_ini to n_end
    def __init__(self, n_ini, n_end, m_ini, alpha_min, alpha_max, alpha_samples):
            
        self.alpha = np.linspace(alpha_min, alpha_max, alpha_samples, endpoint=False)

        # Select shells for initial energy
        # Remember n >= m, and m=n, ..., N, where N it the max wavenumber (space_res_spectral)
        # n defines the shell
        self.nmodes=[]
        self.mmodes=[]
        self.ampls=[]
        self.n_ini = n_ini
        self.n_end = n_end
        self.m_ini = m_ini

        count_modes = 0
        code=""
        
        for n in range(n_ini, n_end+1):
            for m in range(m_ini, n+1):
                self.nmodes.append(n)
                self.mmodes.append(m)
                self.ampls.append(1.0)
                count_modes+=1
                

        self.count_modes = count_modes 

        codes = []
        print()
        print("Mode init params:")
        for a in self.alpha:
            print()
            print("alpha = ", a)
            print("i n m amp")
            code = str(self.count_modes)
            for i in range(self.count_modes):
                code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(a*self.ampls[i])
                print(i, self.nmodes[i], self.mmodes[i], a*self.ampls[i])
            codes.append(code)
        
        self.codes = codes
        print(codes)

    def save_file(self, filename):


        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

class modes_TC2: #list of initial modes
    def __init__(self, n_list, m_list, alpha_min, alpha_max, alpha_samples, back_n_min=0, back_n_max=0, back_ampl=0.1):
            
        self.alpha = np.linspace(alpha_min, alpha_max, alpha_samples, endpoint=False)

        # Select shells for initial energy
        # Remember n >= m, and m=n, ..., N, where N it the max wavenumber (space_res_spectral)
        # n defines the shell
        self.nmodes=n_list
        self.mmodes=m_list
        self.ampls=[]
        self.n_ini = min(n_list)
        self.n_end = max(n_list)
        self.m_ini = min(m_list)

        count_modes = 0
        code=""
        
        for n in n_list:
            self.ampls.append(1.0)
            count_modes+=1
                
        self.count_modes = count_modes 
        list_modes = count_modes

        #add energy on other modes (background energy)
        n_ini = back_n_min
        n_end = back_n_max
        m_ini = 0
        
        if n_ini != 0 and n_end != 0:
            for n in range(n_ini, n_end+1):
                for m in range(m_ini, n+1):
                    if (n,m) in zip(n_list, m_list):
                        continue
                    else:
                        self.nmodes.append(n)
                        self.mmodes.append(m)
                        self.ampls.append(back_ampl)
                        count_modes+=1
                
        self.count_modes = count_modes 

        codes = []
        print()
        print("Mode init params:")
        for a in self.alpha:
            print()
            print("alpha = ", a)
            print("i n m amp")
            code = str(self.count_modes)
            for i in range(self.count_modes):
                if i < list_modes:
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(a*self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], a*self.ampls[i])
                else:
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], self.ampls[i])
            codes.append(code)
        
        self.codes = codes
        print(codes)
        

    def save_file(self, filename):

        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

class modes_TC3: #list of initial modes and list of background modes
    def __init__(self, n_list, m_list, n_list_back, m_list_back, alpha_min, alpha_max, alpha_samples, back_n_min=0, back_n_max=0, back_ampl=0.1):
            
        self.alpha = np.linspace(alpha_min, alpha_max, alpha_samples, endpoint=False)

        # Select shells for initial energy
        # Remember n >= m, and m=n, ..., N, where N it the max wavenumber (space_res_spectral)
        # n defines the shell
        self.nmodes=n_list+n_list_back
        self.mmodes=m_list+m_list_back
        self.ampls=[]
        self.n_ini = min(n_list)
        self.n_end = max(n_list)
        self.m_ini = min(m_list)

        count_modes = 0
        code=""
        
        for n in n_list:
            self.ampls.append(1.0)
            count_modes+=1
                
        self.count_modes = count_modes 
        list_modes = count_modes

        for n in n_list_back:
            self.ampls.append(1.0)
            count_modes+=1

        self.count_modes = count_modes 

        #add energy on other modes (background energy)
        n_ini = back_n_min
        n_end = back_n_max
        m_ini = 0
        
        if n_ini != 0 and n_end != 0:
            for n in range(n_ini, n_end+1):
                for m in range(m_ini, n+1):
                    if (n,m) in zip(self.nmodes, self.mmodes):
                        continue
                    else:
                        self.nmodes.append(n)
                        self.mmodes.append(m)
                        self.ampls.append(back_ampl)
                        count_modes+=1
                
        self.count_modes = count_modes 

        codes = []
        print()
        print("Mode init params:")
        for a in self.alpha:
            print()
            print("alpha = ", a)
            print("i n m amp")
            code = str(self.count_modes)
            for i in range(self.count_modes):
                if i < list_modes:
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(a*self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], a*self.ampls[i])
                else:
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], self.ampls[i])
            codes.append(code)
        
        self.codes = codes
        print(codes)
        

    def save_file(self, filename):

        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

#Read a test case object for post-processing
def load_file(filename):
    f = open(filename, 'rb')
    obj = pickle.load(f)
    f.close()        
    return obj
