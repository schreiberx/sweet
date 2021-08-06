# ---------------------------------------------
# Class to setup spherical modes initialization 
# author: Pedro Peixoto <ppeixoto@usp.br>
# ----------------------------------------
import numpy as np
import pickle

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

def load_file(self, filename):
    f = open(filename, 'rb')
    obj = pickle.load(f)
    f.close()        
    return obj
