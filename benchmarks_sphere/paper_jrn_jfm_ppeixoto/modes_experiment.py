# ---------------------------------------------
# Class to setup spherical modes initialization 
# author: Pedro Peixoto <ppeixoto@usp.br>
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


class modes_TC3_old: ## (legacy - should not be used) - range of modes in shells with added background energy
    def __init__(self, alpha_min, alpha_max, alpha_samples):
            
        self.alpha = np.linspace(alpha_min, alpha_max, alpha_samples, endpoint=False)

        # Select shells for initial energy
        # Remember n >= m, and m=n, ..., N, where N it the max wavenumber (space_res_spectral)
        # n defines the shell
        n_list = [3, 7, 7, 5]
        m_list = [0, 4, 3, 4]
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

        codes = []
        print()
        print("Mode init params:")
        for a in self.alpha:
            print()
            print("alpha = ", a)
            print("i n m amp")
            code = str(self.count_modes)
            for i in range(self.count_modes):
                if i < 2: #first 2 modes get alpha multiplied
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(a*self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], a*self.ampls[i])
                else: #last 2 modes get a small constant
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(0.1)
                    print(i, self.nmodes[i], self.mmodes[i], 0.1)
            codes.append(code)
        
        self.codes = codes
        print(codes)

    def save_file(self, filename):

        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

class modes_TC4_old: # (legacy - should not be used)
    def __init__(self, alpha_min, alpha_max, alpha_samples):
            
        self.alpha = np.linspace(alpha_min, alpha_max, alpha_samples, endpoint=False)

        # Select shells for initial energy
        # Remember n >= m, and m=n, ..., N, where N it the max wavenumber (space_res_spectral)
        # n defines the shell        
        n_list = [5, 3]
        m_list = [4, 1]
        self.nmodes=n_list
        self.mmodes=m_list
        self.ampls=[]
        self.n_ini = min(n_list)
        self.n_end = max(n_list)
        self.m_ini = min(m_list)

        count_modes = 0
        
        for n in n_list:
            self.ampls.append(1.0)
            count_modes+=1
                
        self.count_modes = count_modes 

        #add energy on other modes
        n_ini = 2
        n_end = 4
        m_ini = 0
        
        for n in range(n_ini, n_end+1):
            for m in range(m_ini, n+1):
                if (n,m) in zip(n_list, m_list):
                    continue
                else:
                    self.nmodes.append(n)
                    self.mmodes.append(m)
                    self.ampls.append(0.1)
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
                if i < 2: #first 2 modes get alpha multiplied
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(a*self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], a*self.ampls[i])
                else: #last 2 modes get a small constant
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], 0.1)
            codes.append(code)
        
        self.codes = codes
        print(codes)

    def save_file(self, filename):

        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

class modes_TC5: #list + adds noise (legacy - should not be used)
    def __init__(self, n_list, m_list, alpha_min, alpha_max, alpha_samples):
            
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
                
        #add energy on other modes
        n_ini = 3
        n_end = 5
        m_ini = 0
        
        for n in range(n_ini, n_end+1):
            for m in range(m_ini, n+1):
                if (n,m) in zip(n_list, m_list):
                    continue
                else:
                    self.nmodes.append(n)
                    self.mmodes.append(m)
                    self.ampls.append(0.1)
                    count_modes+=1

        codes = []
        print()
        print("Mode init params:")
        for a in self.alpha:
            print()
            print("alpha = ", a)
            print("i n m amp")
            code = str(count_modes)
            for i in range(count_modes):
                if i < 2: #first 2 modes get alpha multiplied
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(a*self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], a*self.ampls[i])
                else: #last modes get a small constant (noise)
                    code+="_"+str(self.nmodes[i])+"_"+str(self.mmodes[i])+"_"+str(self.ampls[i])
                    print(i, self.nmodes[i], self.mmodes[i], self.ampls[i])
            codes.append(code)
        
        
                
        self.count_modes = count_modes 

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
        self.energy_file    = basedir+"/output_spec_kin_en_t00000000000.00000000.txt"    
        self.enstrophy_file = basedir+"/output_spec_enstrophy_t00000000000.00000000.txt"
        self.scalesmin = {}
        self.scalesmax = {}
        timerescale=1.0/24.00 #Days
        self.energy_file_clean    = basedir+"/output_spec_kin_en_clean_eps"+str(eps)+".pkl"    
        self.enstrophy_file_clean = basedir+"/output_spec_enstrophy_clean_eps"+str(eps)+".pkl"

        #Remove modes with null values

        if os.path.isfile(self.energy_file_clean):
            self.df_energy_clean = pd.read_pickle(self.energy_file_clean)
            maxenergy = self.df_energy_clean['SpectralSum'].max()
            eps_en=maxenergy*eps
            self.scalesmin[0]=eps_en
            self.scalesmax[0]=maxenergy
        else:
            self.df_energy=pd.read_csv(self.energy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
            self.df_energy['timestamp'] = self.df_energy['timestamp'] * timerescale
            maxenergy = self.df_energy['SpectralSum'].max()
            eps_en=maxenergy*eps
            self.scalesmin[0]=eps_en
            self.scalesmax[0]=maxenergy
            self.df_energy_clean = self.df_energy.loc[:, (self.df_energy > eps_en).any(axis=0)]
            self.df_energy_clean = self.df_energy_clean.set_index("timestamp")
            pd.to_pickle(self.df_energy_clean, self.energy_file_clean)

        #print(self.df_energy)
        #df_energy.set_index('timestamp',drop=True,inplace=True)
        #print(self.df_energy_clean)

        if os.path.isfile(self.enstrophy_file_clean):
            self.df_ens_clean = pd.read_pickle(self.enstrophy_file_clean)
            maxens = self.df_ens_clean['SpectralSum'].max()
            eps_ens=maxens*eps
            self.scalesmin[1]=eps_ens
            self.scalesmax[1]=maxens
        else:
            self.df_ens=pd.read_csv(self.enstrophy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
            self.df_ens['timestamp'] = self.df_ens['timestamp'] * timerescale
            maxens = self.df_ens['SpectralSum'].max()
            eps_ens=maxens*eps
            self.scalesmin[1]=eps_ens
            self.scalesmax[1]=maxens
            self.df_ens_clean = self.df_ens.loc[:, (self.df_ens > eps_ens).any(axis=0)]
            self.df_ens_clean = self.df_ens_clean.set_index("timestamp")
            pd.to_pickle(self.df_ens_clean, self.enstrophy_file_clean)
            #print(self.df_ens)
            #df_energy.set_index('timestamp',drop=True,inplace=True)
            #print(self.df_ens_clean)


    def set_out_modes(self, n_list, m_list):

        self.out_modes_name = "out_n"+"-".join([str(i) for i in n_list])+"m"+"-".join([str(i) for i in m_list])
        
        self.df_energy_agg = self.set_out_modes_df(self.df_energy_clean, n_list, m_list)
        
        self.df_ens_agg = self.set_out_modes_df(self.df_ens_clean, n_list, m_list)

        
        self.max_exchange_out_energy = self.df_energy_agg["out_modes"].max()/self.df_energy_agg["init"].iloc[0]
        self.max_exchange_noninit_energy = self.df_energy_agg["non_init"].max()/self.df_energy_agg["init"].iloc[0]

        self.max_exchange_out_ens = self.df_ens_agg["out_modes"].max()/self.df_ens_agg["init"].iloc[0]
        self.max_exchange_noninit_ens = self.df_ens_agg["non_init"].max()/self.df_ens_agg["init"].iloc[0]
        return self.out_modes_name

    def set_out_modes_df(self, df, n_list, m_list):

        tol = 10e-16

        #get modes of initial energy 
        init_modes = []
        non_init_modes = []
        out_modes = []

        ninit_min = 9999
        ninit_max = 0
        for col in df.columns:
            
            if col == "timestamp" or col == "SpectralSum":
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

            mode = (n,m)
            list_modes = list(zip(n_list, m_list))
            #print(list_modes)
            if mode in list_modes:
                out_modes.append(col)
                #print(col, n, m, "out")

        df_init = df[init_modes].sum(axis=1).rename("init") #_n"+str(ninit_min)+"-"+str(ninit_max))
        df_noninit = df[non_init_modes].sum(axis=1).rename("non_init")
        df_out = df[out_modes].sum(axis=1).rename("out_modes")
        #print(df)
        #print(df_out)
        df_new = pd.concat([df_init, df_noninit, df_out], axis=1)
        
        #df_new = df_new.truncate(after=t_limit)
        return df_new

    def set_out_shells(self, nmin, nmax):
        self.out_modes_name = "out_n"+str(nmin)+"-"+str(nmax)
        self.df_energy_agg = self.set_out_shells_df(self.df_energy_clean, nmin, nmax)
        self.df_ens_agg = self.set_out_shells_df(self.df_ens_clean, nmin, nmax)

        
        self.max_exchange_out_energy = self.df_energy_agg["out_modes"].max()/self.df_energy_agg["init"].iloc[0]
        self.max_exchange_noninit_energy = self.df_energy_agg["non_init"].max()/self.df_energy_agg["init"].iloc[0]

        self.max_exchange_out_ens = self.df_ens_agg["out_modes"].max()/self.df_ens_agg["init"].iloc[0]
        self.max_exchange_noninit_ens = self.df_ens_agg["non_init"].max()/self.df_ens_agg["init"].iloc[0]
        
        return self.out_modes_name

    def set_out_shells_df(self, df, nmin, nmax):

        tol = 10e-16

        #get modes of initial energy 
        init_modes = []
        non_init_modes = []
        out_modes = []

        ninit_min = 9999
        ninit_max = 0
        for col in df.columns:
            
            if col == "timestamp" or col == "SpectralSum":
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

        df_init = df[init_modes].sum(axis=1).rename("init") #_n"+str(ninit_min)+"-"+str(ninit_max))
        df_noninit = df[non_init_modes].sum(axis=1).rename("non_init")
        df_out = df[out_modes].sum(axis=1).rename("out_modes")
        
        df_new = pd.concat([df_init, df_noninit, df_out], axis=1)

        return df_new

    def fourier_modes(self, title="", output_filename="out.pdf", do_plot=True, lim_inf=0, lim_sup=999999):

        #Energy
        outmodes = self.df_energy_agg["out_modes"].values
        time = self.df_energy_agg.index.to_numpy()
        #t_limit=-1
        #self.large_periods_energy = self.fourier_largest(outmodes[:t_limit], T=time[t_limit], title=title, output_filename=output_filename, do_plot=do_plot)
        self.fourier_spec = self.fourier_low_freq(outmodes, T=time[-1], title=title, output_filename=output_filename, do_plot=do_plot, lim_inf=lim_inf, lim_sup=lim_sup)
        #Enstrophy
        #outmodes = self.df_energy_agg["out_modes"].values
        #self.large_periods_ens = fourier_largest(outmodes, T=time[-1])

        return self.fourier_spec

    
    def plot(self, title="", output_filename="out.pdf"):

        
        fig, ax = plt.subplots(1, figsize=(6,4))#, sharex=True)
        plt.rc('text', usetex=False)
        
        ax.set_title(title,fontsize=14)

        #for i, ax in enumerate(axs):

        ncol = 1	
        self.df_energy_clean = self.df_energy_clean.reindex((sorted(self.df_energy_clean.columns, reverse=True)), axis=1)
        
        plot_percent = True
        if plot_percent:
            df_percent = self.df_energy_clean.div(self.df_energy_clean['SpectralSum'], axis=0) 
            df_percent = df_percent.drop(columns=['SpectralSum'])
            #print(df_percent)
        else:
            df_percent = self.df_energy_clean

        eps = 0.001
        filter = df_percent.max() > eps
        #print(df_percent)
        df_percent = df_percent.loc[:,filter]
        #print(df_percent)

        color_dict = {}
        if "(5;4) (3;1) (7;3)" in title:
            color_dict = {'(5;4)': 'blue', '(3;1)': 'green', '(7;3)': 'orange', '(9;2)':'red', 'SpectralSum':'gray',
                          '(5;5)': 'blue', '(5;3)': 'green', '(5;1)': 'orange', '(3;2)': 'red'}             
            style_dict = {'(5;4)': '-', '(3;1)': '-', '(7;3)': '-', '(9;2)':'-'}
            linewidths_dict = {'(5;4)': '1', '(3;1)': '1', '(7;3)': '1', '(9;2)':'2'}
            lws = [linewidths_dict.get(x, '0.4') for x in df_percent.columns]

        #self.df_energy_clean.plot( ax=ax, color=[color_dict.get(x, '#333333') for x in self.df_energy_clean.columns])
        df_percent.plot( ax=ax, 
            style=[style_dict.get(x, '--') for x in df_percent.columns], 
            color=[color_dict.get(x, 'gray') for x in df_percent.columns])
        
        for i, l in enumerate(ax.lines):
            plt.setp(l, linewidth=lws[i])
        
        ax.set_yscale("log", nonpositive='clip')
        if plot_percent:
            ax.set(ylabel='$\%$ Energy', xlabel="Time (days)")
            ax.set_ylim([10e-5, 1])
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, 2, '%'))
        else:
            ax.set(ylabel='Energy', xlabel="Time (days)")
        ax.legend(loc='upper left', bbox_to_anchor= (1.0, 1.0), ncol=ncol, fontsize=8)
        ax.set_xscale("linear")
        
        #ylim=[self.scalesmin[0], self.scalesmax[0]]
        #ax.set_ylim(ylim)        

        #fig.subplots_adjust(right=0.7)
        
        print("    ", self.basedir+"/energy_"+output_filename)
        plt.tight_layout()
        plt.savefig(self.basedir+"/energy_"+output_filename, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()

        #enstrophy
        fig, ax = plt.subplots(1, figsize=(10,6))#, sharex=True)
        plt.rc('text', usetex=False)
        
        fig.suptitle(title)

        #for i, ax in enumerate(axs):
        ax.set_xscale("linear")
        ax.set_yscale("log", nonpositive='clip')
        ylim=[self.scalesmin[1], self.scalesmax[1]]
        ax.set_ylim(ylim)

        
        ncol = 1	
    
        self.df_ens_clean.plot( ax=ax)
        ax.set(ylabel='Enstrophy', xlabel="Time (days)")
        #axs[0].get_legend().remove()
        ax.legend(loc='upper left', bbox_to_anchor= (1.0, 1.0), ncol=ncol, fontsize="small")

        fig.subplots_adjust(right=0.7)
        
        print("    ", self.basedir+"/enstrophy_"+output_filename)
        #plt.tight_layout()
        plt.savefig(self.basedir+"/enstrophy_"+output_filename, transparent=True, dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()

    def plot_out(self, title="", output_filename="out.png"):

        fontsize=18
        figsize=(10, 10)

        fig, axs = plt.subplots(2, figsize=(10,10))#, sharex=True)
        plt.rc('text', usetex=False)
        #title="SPH Mode Nonlinear Interaction\n"+title
        
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

        
        ncol = 1	
        
        self.df_energy_agg.plot( ax=axs[0])
        axs[0].set(ylabel='Energy', xlabel="Time (days)")
        axs[0].legend(loc='center left',  bbox_to_anchor= (1.0, -0.1), ncol=ncol, fontsize="small")
        
        ncol=1
    
        self.df_ens_agg.plot( ax=axs[1])
        axs[1].set(ylabel='Enstrophy', xlabel="Time (days)")
        axs[1].get_legend().remove()
        #axs[1].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

        fig.subplots_adjust(right=0.7)
        plt.tight_layout()

        print("    ", self.basedir+"/"+output_filename)
        #plt.show()
        plt.savefig(self.basedir+"/"+output_filename, transparent=True, dpi=300) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()


    def fourier_largest(self, array, T, title="", output_filename="out.pdf", do_plot=True):
        #plt.plot(array)
        #plt.show()
        if np.max(np.abs(array)) < 0.000001:
            return [0]
        yf = np.fft.rfft(array)
        #cut out zero freq
        yf = np.abs(yf[1:])
        yf = np.abs(yf)**2 #power spectrum
        #plt.plot(yf)
        #plt.show()
        n = len(yf)

        xf = T/np.linspace(1, T, n)

        nlargest = 10
        ilargest = np.argpartition(yf, -nlargest)[-nlargest:]
        large_periods = np.take(xf, ilargest)
        large_spectrum = np.take(yf, ilargest)

        #Ordered largest elements of spectrum
        large_sort = large_spectrum.argsort()
        large_periods = large_periods[large_sort[::-1]]
        large_spectrum = large_spectrum[large_sort[::-1]]
        print(large_periods)
        print(large_spectrum)

        #Filter modes of very large period
        #large_filter = large_periods<T/2
        #large_periods = large_periods[large_filter]
        #large_spectrum = large_spectrum[large_filter]
        #print(large_periods)
        #print(large_spectrum)

                
        #print(xf)
        #print(yf)
        if do_plot:
            fig, ax = plt.subplots()
            #ax.set_xscale('log')
            fig.suptitle(title)
            #ax.plot(xf, 1.0/n * yf)
            ax.plot(xf, yf)
            ax.set(ylabel='Power Spectrum', xlabel="Periodicity (days)")
            #ax.set(ylabel='Normalized Mode Amplitude', xlabel="Periodicity (days)", title=title)
            plt.xscale('log', base=10)
            plt.yscale('log', base=10)
            #ax.plot(yf)
            plt.tight_layout()
            print("    ", self.basedir+"/"+output_filename)
            plt.savefig(self.basedir+"/"+output_filename, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)
            plt.close()

        return large_periods

    def fourier_low_freq(self, array, T, title="", output_filename="out.pdf", do_plot=True, lim_inf=0, lim_sup=99999):
        #plt.plot(array)
        #plt.show()
        if np.max(np.abs(array)) < 0.000001:
            return 0
        yf = np.fft.rfft(array)
        #cut out zero freq
        yf = np.abs(yf[1:])
        yf = np.abs(yf) #power spectrum
        #plt.plot(yf)
        #plt.show()
        n = len(yf)

        xf = T/np.linspace(1, T, n)

        #Filter modes 
        filter1 = xf > lim_inf
        filter2 = xf < lim_sup
        filter = filter1 & filter2
        
        #filter =  lim_inf < xf and xf < lim_sup
        periods = xf[filter]
        spectrum = yf[filter]
       
        power_spec_filtred = spectrum.sum()

        #print(periods)
        #print(spectrum)
        #print(power_spec_filtred, np.any(filter))

        #print(xf)
        #print(yf)
        if do_plot:
            fig, ax = plt.subplots(figsize=(5,5))
            
            ax.set_title(title, fontsize=14)
            #ax.plot(xf, 1.0/n * yf)
            ax.plot(xf[2:], yf[2:]**2, color='red')
            ax.set(ylabel='Power Spectrum', xlabel="Periodicity (days)")
            #ax.set(ylabel='Normalized Mode Amplitude', xlabel="Periodicity (days)", title=title)
            plt.xscale('log', base=10)
            plt.yscale('log', base=10)
            ax.set_ylim([10e-5, 10e3])
            #ax.plot(yf)
            plt.tight_layout()
            print("    ", self.basedir+"/"+output_filename)
            plt.savefig(self.basedir+"/"+output_filename, transparent=True, dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)
            plt.close()

        return power_spec_filtred
