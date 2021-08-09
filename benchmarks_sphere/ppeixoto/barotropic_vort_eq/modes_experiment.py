# ---------------------------------------------
# Class to setup spherical modes initialization 
# author: Pedro Peixoto <ppeixoto@usp.br>
# ----------------------------------------
import numpy as np
import pickle
from numpy.lib.function_base import append
import pandas as pd
import re

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D

from mule.postprocessing.JobData import *

class modes:
    def __init__(self, n_ini, n_end, m_ini, alpha_min, alpha_max, alpha_samples):
            
        self.alpha = np.linspace(alpha_min, alpha_max, alpha_samples, endpoint=True)

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

def load_file(filename):
    f = open(filename, 'rb')
    obj = pickle.load(f)
    f.close()        
    return obj

class evol:
    def __init__(self, basedir=".", eps=0.0001):

        self.basedir = basedir
        self.energy_file    = basedir+"/output_spec_energy_t00000000000.00000000.txt"    
        self.enstrophy_file = basedir+"/output_spec_enstrophy_t00000000000.00000000.txt"
        self.scalesmin = {}
        self.scalesmax = {}
        timerescale=1.0

        #Remove modes with null values

        self.df_energy=pd.read_csv(self.energy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
        #print(self.df_energy)
        self.df_ens=pd.read_csv(self.enstrophy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
        #print(self.df_ens)

        self.df_energy['timestamp'] = self.df_energy['timestamp'] * timerescale
        maxenergy = self.df_energy['TotalSum'].max()
        eps_en=maxenergy*eps
        self.scalesmin[0]=eps_en
        self.scalesmax[0]=maxenergy
        self.df_energy_clean = self.df_energy.loc[:, (self.df_energy > eps_en).any(axis=0)]
        self.df_energy_clean = self.df_energy_clean.set_index("timestamp")
        #print(self.df_energy)
        #df_energy.set_index('timestamp',drop=True,inplace=True)
        print(self.df_energy_clean)

        self.df_ens['timestamp'] = self.df_ens['timestamp'] * timerescale
        maxens = self.df_ens['TotalSum'].max()
        eps_ens=maxens*eps
        self.scalesmin[1]=eps_ens
        self.scalesmax[1]=maxens
        self.df_ens_clean = self.df_ens.loc[:, (self.df_ens > eps_ens).any(axis=0)]
        self.df_ens_clean = self.df_ens_clean.set_index("timestamp")
        #print(self.df_ens)
        #df_energy.set_index('timestamp',drop=True,inplace=True)
        print(self.df_ens_clean)


    def set_out_shells(self, nmin, nmax):
        self.df_energy_agg = self.set_out_shells_df(self.df_energy_clean, nmin, nmax)
        self.df_ens_agg = self.set_out_shells_df(self.df_ens_clean, nmin, nmax)

    def set_out_shells_df(self, df, nmin, nmax):

        tol = 10e-16

        #get modes of initial energy 
        init_modes = []
        non_init_modes = []
        out_modes = []

        ninit_min = 9999
        ninit_max = 0
        for col in df.columns:
            
            if col == "timestamp" or col == "TotalSum":
                continue
            a = [int(s) for s in re.findall(r'\b\d+\b', col)]
            n = a[0]
            m = a[1]
            #print(df[col].iloc[0])
            if df[col].iloc[0] > tol:
                init_modes.append(col)
                #print(col, n, m, "init")
                if n > ninit_max:
                    ninit_max = n
                if n < ninit_min:
                    ninit_min = n
            else:
                non_init_modes.append(col)
                #print(col, n, m, "non-init")

            if n >= nmin and n <= nmax:
                out_modes.append(col)
                #print(col, n, m, "out")
            
        df_init = df[init_modes].sum(axis=1).rename("init_n"+str(ninit_min)+"-"+str(ninit_max))
        df_noninit = df[non_init_modes].sum(axis=1).rename("non_init")
        df_out = df[out_modes].sum(axis=1).rename("out_n"+str(nmin)+"-"+str(nmax))
        
        df_new = pd.concat([df_init, df_noninit, df_out], axis=1)

        return df_new

    def plot(self, title="", output_filename="out.pdf"):

        
        fontsize=18
        figsize=(10, 10)

        fig, axs = plt.subplots(2, figsize=(10,10))#, sharex=True)
        plt.rc('text', usetex=False)
        title="SPH Mode Nonlinear Interaction\n"+title
        
        fig.suptitle(title)

        
        for i, ax in enumerate(axs):
            ax.set_xscale("linear")
            ax.set_yscale("log", nonpositive='clip')
            ylim=[self.scalesmin[i], self.scalesmax[i]]
            ax.set_ylim(ylim)

        #for ax in axs.flat:
        #    ax.label_outer()

        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

        markers = []
        for m in Line2D.markers:
            try:
                if len(m) == 1 and m != ' ' and m != '':
                    markers.append(m)
            except TypeError:
                pass

        linestyles = ['-', '--', ':', '-.']

        
        ncol = 2	
        
        self.df_energy_clean.plot( ax=axs[0])
        axs[0].set(ylabel='Energy', xlabel="Time (hours)")
        axs[0].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)
        
        ncol=2
    
        self.df_ens_clean.plot( ax=axs[1])
        axs[1].set(ylabel='Enstrophy', xlabel="Time (hours)")
        axs[1].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

        fig.subplots_adjust(right=0.7)
        
        print(self.basedir+"/"+output_filename)
        plt.show()
        plt.savefig(self.basedir+"/"+output_filename, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()

    def plot_shells(self, title="", output_filename="out.pdf"):

        fontsize=18
        figsize=(10, 10)

        fig, axs = plt.subplots(2, figsize=(10,10))#, sharex=True)
        plt.rc('text', usetex=False)
        title="SPH Mode Nonlinear Interaction\n"+title
        
        fig.suptitle(title)

        for i, ax in enumerate(axs):
            ax.set_xscale("linear")
            ax.set_yscale("log", nonpositive='clip')
            ylim=[self.scalesmin[i], self.scalesmax[i]]
            ax.set_ylim(ylim)

        #for ax in axs.flat:
        #    ax.label_outer()

        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

        markers = []
        for m in Line2D.markers:
            try:
                if len(m) == 1 and m != ' ' and m != '':
                    markers.append(m)
            except TypeError:
                pass

        linestyles = ['-', '--', ':', '-.']

        
        ncol = 2	
        
        self.df_energy_agg.plot( ax=axs[0])
        axs[0].set(ylabel='Energy', xlabel="Time (hours)")
        axs[0].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)
        
        ncol=2
    
        self.df_ens_agg.plot( ax=axs[1])
        axs[1].set(ylabel='Enstrophy', xlabel="Time (hours)")
        axs[1].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

        fig.subplots_adjust(right=0.7)
        
        print(self.basedir+"/"+output_filename)
        plt.show()
        plt.savefig(self.basedir+"/"+output_filename, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()

