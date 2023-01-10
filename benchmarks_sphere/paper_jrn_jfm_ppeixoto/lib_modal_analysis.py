# ---------------------------------------------
# Class to setup spherical modes post-processing
# author: Pedro Peixoto <ppeixoto@usp.br>
#  Oct 2021
# ----------------------------------------
import numpy as np
import pickle
from numpy.lib.function_base import append
import pandas as pd
import re
import os
import os.path
from scipy import stats
from scipy import signal

import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import matplotlib.ticker as mtick


from mule.postprocessing.JobData import *

class plot_fmt:
    def __init__(self, color_dict={}, style_dict={}, linewidths_dict={}, valid=False):
        self.color_dict = color_dict
        self.style_dict = style_dict
        self.linewidths_dict = linewidths_dict
        if len(color_dict) > 0:
            self.valid = True


class mode_evol:
    def __init__(self, basedir=".", eps=0.0001):

        self.basedir = basedir
        #self.energy_file    = basedir+"/output_spec_kin_en_t00000000000.00000000.txt"
        #self.enstrophy_file = basedir+"/output_spec_enstrophy_t00000000000.00000000.txt"
        self.energy_file    = basedir+"/output_spec_ampl_kin_en.txt"    
        self.enstrophy_file = basedir+"/output_spec_ampl_enstrophy.txt"
        self.energy_phase_file    = basedir+"/output_spec_arg_kin_en.txt"    
        self.enstrophy_phase_file = basedir+"/output_spec_arg_enstrophy.txt"
        self.scalesmin = {}
        self.scalesmax = {}
        timerescale=1.0/24.00 #Days
        self.energy_file_clean    = basedir+"/output_spec_kin_en_clean_eps"+str(eps)+".pkl"    
        self.enstrophy_file_clean = basedir+"/output_spec_enstrophy_clean_eps"+str(eps)+".pkl"
        self.energy_phase_file_clean    = basedir+"/output_spec_arg_kin_en_clean_eps"+str(eps)+".pkl"    
        self.enstrophy_phase_file_clean = basedir+"/output_spec_arg_enstrophy_clean_eps"+str(eps)+".pkl"

        iphase = True #Analyse phase

        #Remove modes with null values
        if os.path.isfile(self.energy_file_clean):
            self.df_energy_clean = pd.read_pickle(self.energy_file_clean)
            maxenergy = self.df_energy_clean['SpectralSum'].max()
            eps_en = maxenergy*eps
            self.scalesmin[0] = eps_en
            self.scalesmax[0] = maxenergy
        else:
            self.df_energy = pd.read_csv(self.energy_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
            self.df_energy['timestamp'] = self.df_energy['timestamp'] * timerescale
            maxenergy = self.df_energy['SpectralSum'].max()
            eps_en = maxenergy*eps
            self.scalesmin[0] = eps_en
            self.scalesmax[0] = maxenergy
            self.df_energy_clean = self.df_energy.loc[:, (self.df_energy > eps_en).any(axis=0)]
            self.df_energy_clean = self.df_energy_clean.set_index("timestamp")
            pd.to_pickle(self.df_energy_clean, self.energy_file_clean)

        self.df_energy_clean = self.df_energy_clean.reindex((sorted(self.df_energy_clean.columns, reverse=True)), axis=1)

        if iphase and os.path.isfile(self.energy_phase_file_clean):
            self.df_energy_phase_clean = pd.read_pickle(self.energy_phase_file_clean)
        else:            
            if iphase:
                self.df_energy_phase = pd.read_csv(self.energy_phase_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
                self.df_energy_phase = self.df_energy_phase.iloc[:, :-1] #remove garbage column that come out of sweet
                
                self.df_energy_phase['timestamp'] = self.df_energy_phase['timestamp'] * timerescale
                    
                self.df_energy_phase_clean = self.df_energy_phase.loc[:, (self.df_energy > eps_en).any(axis=0)]
                self.df_energy_phase_clean = self.df_energy_phase_clean.set_index("timestamp")
                pd.to_pickle(self.df_energy_phase_clean, self.energy_phase_file_clean)
            

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

        self.df_ens_clean = self.df_ens_clean.reindex((sorted(self.df_ens_clean.columns, reverse=True)), axis=1)

        if iphase and os.path.isfile(self.enstrophy_phase_file_clean):
            self.df_enstrophy_phase_clean = pd.read_pickle(self.enstrophy_phase_file_clean)
        else:
            if iphase:
                self.df_enstrophy_phase = pd.read_csv(self.enstrophy_phase_file, sep='\t', skipinitialspace=True, skiprows=1, header=3, engine="python")
                self.df_enstrophy_phase = self.df_enstrophy_phase.iloc[:, :-1] #remove garbage column that come out of sweet
                
                self.df_enstrophy_phase['timestamp'] = self.df_enstrophy_phase['timestamp'] * timerescale
                    
                self.df_enstrophy_phase_clean = self.df_enstrophy_phase.loc[:, (self.df_ens > eps_ens).any(axis=0)]
                self.df_enstrophy_phase_clean = self.df_enstrophy_phase_clean.set_index("timestamp")
                pd.to_pickle(self.df_enstrophy_phase_clean, self.enstrophy_phase_file_clean)

        
        


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

    def phase(self):

        #print(self.df_energy_phase_clean)
        df = self.df_energy_phase_clean.copy()
        t = df.index.values
        dt = t[1]-t[0]
        #print(dt)
        for column in df:
            x = np.unwrap(df[column].values)
            dif = (x[1:]-x[0:-1])/dt
            dif = np.insert(dif, 0, 0.0, axis=0)
            df[column] = dif  #differential d/dt
            #df[column] = x #show just unwrapped phases
        self.df_energy_phase_dif = df.iloc[4:] #Remove first hours (mode adjustment)
        #print(df)
        #self.df_energy_phase_dif['(9;2)']=self.df_energy_phase_dif['(9;2)'].rolling(72).mean()
        #df_tmp = self.df_energy_phase_dif['(9;2)']
        #self.df_energy_phase_dif['(9;2)'] = df_tmp[(np.abs(stats.zscore(df_tmp)) < 2)]
        self.df_dif_mean = df.mean(axis=0)
        return self.df_dif_mean

    def fourier_low_freq(self, array, T, lim_inf=0, lim_sup=99999):
        
        fs = T/(len(array)-1) #1 hour frequency
        
        if np.max(np.abs(array)) < 0.000001:
            return 0, [0]
        yf = np.fft.rfft(array)
        
        #cut out zero freq
        yf = np.abs(yf[1:])
        yf = np.abs(yf) #power spectrum
        #plt.plot(yf)
        #plt.show()
        n = len(yf)
        #https://en.wikipedia.org/wiki/Spectral_density
        dt = 1 #in hours #1/24 hours -> per day
        pds = (dt*dt*yf**2) /T #all is in hours
        xf = T/np.linspace(1, T, n)
        #print(xf)
        #convert to days for filter
        xf = xf / 24
        #print(xf[1:5])
        
        #print(xf)
        #Filter modes 
        filter1 = xf > lim_inf
        filter2 = xf < lim_sup
        filter = filter1 & filter2
        
        
        #filter =  lim_inf < xf and xf < lim_sup
        #periods = xf[filter]
        spectrum = yf[filter]
       
        power_spec_filtred = dt*dt*np.sqrt((spectrum**2).sum()) /T

        #save spectrum
        #self.spec = yf

        return power_spec_filtred, pds

    #-----------------------------------------------
    #  Ploting functions
    #------------------------------------------------



    def plot_modes(self, df, var, title="", output_filename="out.pdf", pltfmt=None):

        
        fig, ax = plt.subplots(1, figsize=(6,4))#, sharex=True)
        plt.rc('text', usetex=False)
        
        ax.set_title(title,fontsize=14)

        #for i, ax in enumerate(axs):        
        
        plot_percent = True
        if plot_percent:
            df_percent = df.div(df['SpectralSum'], axis=0) 
            df_percent = df_percent.drop(columns=['SpectralSum'])
            #print(df_percent)
        else:
            df_percent = df

        if var == "Energy":
            eps = 0.001
        else:
            eps = 0.01

        filter = df_percent.max() > eps
        #print(df_percent)
        df_percent = df_percent.loc[:,filter]
        #print(df_percent)

        if pltfmt is not None:
            color_dict = pltfmt.color_dict 
            style_dict = pltfmt.style_dict 
            linewidths_dict = pltfmt.linewidths_dict 
        else:
            color_dict = {}
            style_dict = {}
            linewidths_dict = {}
        lws = [linewidths_dict.get(x, '0.4') for x in df_percent.columns]

        if len(color_dict)>0:
            df_percent.plot( ax=ax, 
                style=[style_dict.get(x, ':') for x in df_percent.columns], 
                color=[color_dict.get(x, 'gray') for x in df_percent.columns])
            for i, l in enumerate(ax.lines):
                plt.setp(l, linewidth=lws[i])
        else:
            df_percent.plot( ax=ax)

        ax.set_yscale("log", nonpositive='clip')
        if plot_percent:
            ax.set(ylabel='$\%$ '+var, xlabel="Time (days)")
            ax.set_ylim([10e-5, 1])
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, 2, '%'))
        else:
            ax.set(ylabel=var, xlabel="Time (days)")
        ax.legend(loc='upper left', bbox_to_anchor= (1.0, 1.0), fontsize=8)
        ax.set_xscale("linear")
        
        #ylim=[self.scalesmin[0], self.scalesmax[0]]
        #ax.set_ylim(ylim)        

        #fig.subplots_adjust(right=0.7)
        
        print("    ", self.basedir+"/"+output_filename)
        plt.tight_layout()
        plt.savefig(self.basedir+"/"+output_filename, transparent=True, dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()

    def plot_total(self, df, var, title="", output_filename="out.pdf", pltfmt=None):

        
        fig, ax = plt.subplots(1, figsize=(6,4))#, sharex=True)
        plt.rc('text', usetex=False)
        
        ax.set_title(title,fontsize=14, y=1.08)
        
        init_en = df['SpectralSum'].iloc[0]
        df_percent = (df['SpectralSum']-init_en)/init_en
        
        df_percent.plot( ax=ax)

        #ax.set_yscale("log", nonpositive='clip')
        
        ax.set(ylabel='$\% $ Var in '+var, xlabel="Time (days)")
        #ax.set_ylim([10e-5, 1])
        #ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, 2, '%'))
        
        #ax.set(ylabel=var, xlabel="Time (days)")
        #ax.legend(loc='upper left', bbox_to_anchor= (1.0, 1.0), fontsize=8)
        ax.set_xscale("linear")
        
        #ylim=[self.scalesmin[0], self.scalesmax[0]]
        #ax.set_ylim(ylim)        

        #fig.subplots_adjust(right=0.7)
        
        print("    ", self.basedir+"/"+output_filename)
        plt.tight_layout()
        plt.savefig(self.basedir+"/"+output_filename, transparent=True, dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)

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

       
        self.df_energy_agg.plot( ax=axs[0], color = ['b', 'g', 'r'])
        axs[0].set(ylabel='Energy', xlabel="Time (days)")
        axs[0].legend(loc='center left',  bbox_to_anchor= (1.0, -0.1),  fontsize="small")
        
    
        self.df_ens_agg.plot( ax=axs[1], color = ['b', 'g', 'r'])
        axs[1].set(ylabel='Enstrophy', xlabel="Time (days)")
        axs[1].get_legend().remove()
        #axs[1].legend(loc='center left', bbox_to_anchor= (1.01, 0.5), ncol=ncol)

        fig.subplots_adjust(right=0.7)
        plt.tight_layout()

        print("    ", self.basedir+"/"+output_filename)
        #plt.show()
        plt.savefig(self.basedir+"/"+output_filename, transparent=True, dpi=300) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()

    def plot_phase(self, mode_list, title, output_filename="out.pdf"):
        
        fig, ax = plt.subplots(1, figsize=(8,6))#, sharex=True)
        plt.rc('text', usetex=False)
        
        ax.set_title(title,fontsize=14)
        
        cols = set(mode_list)-set(self.df_energy_phase_dif.columns)
        cols = list(set(mode_list) - cols)
        
        df=self.df_energy_phase_dif[list(cols)]
        if df.empty:
            return

        df.plot(kind='line', ax=ax)
        
        ax.set(ylabel='Deriv Phase (rad/s)', xlabel="Time (days)")
        ax.set_ylim([-2, 3])
        print("    ", self.basedir+"/energy_"+output_filename)
        plt.tight_layout()
        plt.savefig(self.basedir+"/energy_"+output_filename, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)

        plt.close()

    def fourier_plot(self, yf, T=365, title="", output_filename="out.pdf"):
        
        if len(yf) < 2:
            return
        n = len(yf)
        xf = T/np.linspace(1, T, n)
        xf = xf / 24 #(convert to hours)
        
        fig, ax = plt.subplots(figsize=(5,5))
        
        ax.set_title(title, fontsize=14)
        #ax.plot(xf, 1.0/n * yf)
        ax.plot(xf, yf, color='red')
        ax.set(ylabel='Power Spectrum Density', xlabel="Periodicity (days)")
        #ax.set(ylabel='Normalized Mode Amplitude', xlabel="Periodicity (days)", title=title)
        plt.xscale('log', base=10)
        plt.yscale('log', base=10)
        ax.set_ylim([10e-7, 10e2])
        ax.set_xlim([1, 370])
        #ax.plot(yf)
        plt.tight_layout()
        print("    ", self.basedir+"/"+output_filename)
        plt.savefig(self.basedir+"/"+output_filename, transparent=True, dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)
        plt.close()

        return 



def plot_uv(df, title, filename_final):

	#Plot velocities
	
	plt.figure(figsize=(10,6), tight_layout=True)
	plt.plot(df, '-', linewidth=2)

	plt.xlabel(r" $\alpha$")
	plt.ylabel(r" Max Velocity $m/s$")
	
	plt.title(title)
	plt.legend(title_fontsize = 13, labels=['zonal (u)' , 'meridional (v)'])

	print("Vel file:", filename_final)
	plt.savefig(filename_final, transparent=True) #, bbox_inches='tight') #, pad_inches=0.02)
	plt.close()


def plot_v_alpha(x, y, z, var, title, spec_string="", filename_final="out.png"):

    fig, ax = plt.subplots(figsize=(5,5)) #, tight_layout=True)
    
    if var == "Energy":
        ax.plot(x, y,  linewidth=1, color = 'blue', linestyle = '-', label=var) #, '--', ':'])
        ax.set_ylabel(r" $\epsilon(\alpha)$", color='blue', fontsize=16)
    else:
        ax.plot(x, y,  linewidth=1, color = 'green', linestyle = '-', label=var) #, '--', ':'])
        ax.set_ylabel(r" $\epsilon(\alpha)$", color='green', fontsize=16)
    ax.set_xlabel(r" $\alpha$", fontsize=16)

    ax.yaxis.set_major_formatter(mtick.PercentFormatter())

    if len(spec_string)>1:
        ax2=ax.twinx()
        ax2.plot(x, z,  linewidth=1, color = 'red', linestyle = ':', label='Low Freq Spec\n'+spec_string) #, '--', ':'])
        ax2.set_ylabel(r'$\sqrt{PSD}$', color='red')

    
    plt.title(title, fontsize=14)

    ax.legend(loc=2, fontsize=12)
    ax2.legend(loc=4, fontsize=11)

    print("Output file:", filename_final)
    plt.tight_layout()
    plt.savefig(filename_final, transparent=True,  dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)
    plt.close()


def plot_omega_v_alpha(x, y, z, w, var, title, spec_string="", filename_final="out.png"):

    fig, ax = plt.subplots(figsize=(5,5)) #, tight_layout=True)
        
    ax.plot(x, y,  linewidth=1, color = 'blue', linestyle = '-', label=var+"-Main") #, '--', ':'])
    ax.set_ylabel(r" $\Omega_N$", fontsize=16)
    
    ax.plot(x, z,  linewidth=1, color = 'green', linestyle = '-', label=var+"-Out") #, '--', ':'])
    ax.set_xlabel(r" $\alpha$", fontsize=16)

    if len(spec_string)>1:
        ax2=ax.twinx()
        ax2.plot(x, w,  linewidth=1, color = 'red', linestyle = ':', label='Low Freq Spec\n'+spec_string) #, '--', ':'])
        ax2.set_ylabel(r'$\sqrt{PSD}$', color='red')
    
    plt.title(title, fontsize=14)

    ax.legend(fontsize=12, bbox_to_anchor=(0.5, 0.8))
    ax2.legend(loc=4, fontsize=11) #, bbox_to_anchor=(1.1, 1.05))

    print("Output file:", filename_final)
    plt.tight_layout()
    plt.savefig(filename_final, transparent=True,  dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)
    plt.close()

def plot_panel(df, title, spec_string="", filename_final="out.png"):

    fig, axs = plt.subplots(4, figsize=(8,8), tight_layout=True)
    axs[0].set_title(title, fontsize=14)
    

    axs[0].plot(df.index, df['Exch_energy']*100,  linewidth=1, color = 'blue', linestyle = '-', label="Energy") #, '--', ':'])
    axs[0].set_ylabel(r"$\epsilon$ (Energy) ", fontsize=14)
    #axs[0].legend(fontsize=12)
    axs[0].yaxis.set_major_formatter(mtick.PercentFormatter())
    if "small" in title:
        axs[0].text(30.0, 1.2, "(a)")
    else:
        axs[0].text(30.0, 2.0, "(a)")

    axs[1].plot(df.index, df['Exch_ens']*100,  linewidth=1, color = 'blue', linestyle = '-', label="Enstrophy") #, '--', ':'])
    axs[1].set_ylabel(r"$\epsilon$ (Enstrophy)", fontsize=13)
    axs[1].yaxis.set_major_formatter(mtick.PercentFormatter())
    #axs[1].legend(fontsize=12, bbox_to_anchor=(0.5, 0.8))
    if "small" in title:
        axs[1].text(30.0, 5.0, "(b)")
    else:
        axs[1].text(30.0, 8.0, "(b)")

    #if len(spec_string)>1:
    #    ax2=ax.twinx()
    #    ax2.plot(x, w,  linewidth=1, color = 'red', linestyle = ':', label='Low Freq Spec\n'+spec_string) #, '--', ':'])
    #    ax2.set_ylabel(r'$\sqrt{PSD}$', color='red')
    
    axs[2].plot(df.index, df['SpecEnerg'],  linewidth=1, color = 'blue', linestyle = '-', label='Low Freq Spec\n'+spec_string) #, '--', ':'])
    #axs[2].set_ylabel("Power Spectrum", fontsize=13)
    axs[2].set_ylabel(r'$\sqrt{PSD}$', fontsize=13)
    axs[2].legend(fontsize=12, loc=2)
    if "small" in title:
        axs[2].text(30.0, 0.028, "(c)")
    else:
        axs[2].text(30.0, 0.06, "(c)")
    axs[2].yaxis.set_major_locator(plt.MaxNLocator(4))
    
    y1 = np.abs(df['Omega_Main'])/3600
    y2 = np.abs(df['Omega_Out'])/3600
    axs[3].plot(df.index, y1,  linewidth=1, color = 'blue', linestyle = '--', label=r" $|\Omega_{12}^3|$") #, '--', ':'])
    axs[3].plot(df.index, y2,  linewidth=1, color = 'orange', linestyle = '-', label=r" $|\Omega_{14}^2|$") #, '--', ':'])
    axs[3].set_ylabel(r" Preces. Freq. ($Hz$)", fontsize=13)
    axs[3].set_xlabel(r" $\alpha$", fontsize=16)
    axs[3].legend(fontsize=12)
    axs[3].yaxis.set_major_locator(plt.MaxNLocator(4))
    #axs[3].yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.1e'))
    #axs[3].yaxis.set_ticks([0.0, 0.00005, 0.00010, 0.00015])
    formatter = mtick.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    axs[3].yaxis.set_major_formatter(formatter) 
    if "small" in title:
        axs[3].text(30.0, 0.00013, "(d)")    
    else:
        axs[3].text(30.0, 0.00016, "(d)")

    #plt.xscale('log', base=10)
    #plt.yscale('log', base=10)
    #ax.set_ylim([10e-7, 10e2])
    #ax.set_xlim([1, 370])
    #ax2.legend(loc=4, fontsize=11) #, bbox_to_anchor=(1.1, 1.05))

    print("Output file:", filename_final)
    plt.tight_layout()
    plt.savefig(filename_final, transparent=True,  dpi=600) #, bbox_inches='tight') #, pad_inches=0.02)
    plt.close()
