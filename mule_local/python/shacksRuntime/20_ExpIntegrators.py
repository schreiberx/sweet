from mule.JobCompileOptions import *


class ExpIntegrators:
    def __init__(self):

        #
        # ExpIntegrators method:
        #
        # '': default
        # 'file': load coefficients from file (terry, cauchy & butcher coefficients)
        # 'direct': Use direct solution, if available
        # 'terry': deprecated T-REXI method
        # 'butcher': deprecated Butcher-REXI method
        #
        self.rexi_method = ''
    
        # Generic REXI parameters
        self.rexi_sphere_preallocation = 0
    
        # List of REXI Coefficients
        self.rexi_files_coefficients = []
    
        # Parameters for T-REXI method
        self.rexi_terry_m = 0
        self.rexi_terry_h = 0.15
        self.rexi_terry_reduce_to_half = 0
    
        # Parameters for Cauchy Contour integral method
        self.rexi_ci_primitive = 'circle'
        self.rexi_ci_n = 128
        self.rexi_ci_max_real = None
        self.rexi_ci_max_imag = None
        self.rexi_ci_sx = None
        self.rexi_ci_sy = None
        self.rexi_ci_mu = None

        # Parameters for direct exp. solver
        self.exp_direct_precompute_phin = 0


    def load_from_dict(self, d):

        if 'rexi_method' in d:
            self.rexi_method = d['rexi_method']

        if 'sphere_preallocation' in d:
            self.rexi_sphere_preallocation = d['sphere_preallocation']

        if 'terry_m' in d:
            self.rexi_terry_m = d['terry_m']

        if 'terry_h' in d:
            self.rexi_terry_h = d['terry_h']


        if 'ci_n' in d:
            self.rexi_ci_n = int(d['ci_n'])

        if 'ci_max_real' in d:
            self.rexi_ci_max_real = float(d['ci_max_real'])

        if 'ci_max_imag' in d:
            self.rexi_ci_max_imag = float(d['ci_max_imag'])

        if 'ci_sx' in d:
            self.rexi_ci_sx = float(d['ci_sx'])

        if 'ci_sy' in d:
            self.rexi_ci_sy = float(d['ci_sy'])

        if 'ci_mu' in d:
            self.rexi_ci_mu = float(d['ci_mu'])

        if 'ci_primitive' in d:
            self.rexi_ci_primitive = float(d['ci_primitive'])

        if 'exp_direct_precompute_phin' in d:
            self.exp_direct_precompute_phin = d['exp_direct_precompute_phin']



    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        idstr = ''

        if not 'runtime.rexi' in filter_list:
            if self.rexi_method != '' and self.rexi_method != None:
                if self.rexi_method == 'direct':
                    idstr += '_REXIDIRECT'
                else:
                    if self.rexi_method == "file":
                        if len(self.rexi_files_coefficients) != 0:
                            idstr += '_FREXI'
                            if not 'runtime.rexi_params' in filter_list:

                                for i in range(len(self.rexi_files_coefficients)):
                                    r = self.rexi_files_coefficients[i]
                                    if r.unique_id_string == None:
                                        raise Exception("Unique ID String missing for REXI coefficients")

                                    if not 'runtime.rexi_params_phi'+str(i) in filter_list:
                                        idstr += "_"+r.unique_id_string.replace("REXI", "")


                    elif self.rexi_method == "terry":
                        idstr += '_TREXI'
                        if not 'runtime.rexi_params' in filter_list:
                            idstr += '_m'+str(self.rexi_terry_m).zfill(8)
                            idstr += '_h'+str(self.rexi_terry_h)
                            idstr += '_r'+str(self.rexi_terry_reduce_to_half)

                    elif self.rexi_method == "ci":
                        idstr += '_CIREXI'

                        if not 'runtime.rexi_params' in filter_list:
                            idstr += '_n'+str(self.rexi_ci_n).zfill(8)
                            if self.rexi_ci_max_real != None:
                                if self.rexi_ci_max_real != None:
                                    idstr += '_mr'+str(float(self.rexi_ci_max_real))
                                if self.rexi_ci_max_imag != None:
                                    idstr += '_mi'+str(float(self.rexi_ci_max_imag))
                            else:
                                if self.rexi_ci_sx != None:
                                    idstr += '_sx'+str(float(self.rexi_ci_sx))
                                if self.rexi_ci_sy != None:
                                    idstr += '_sy'+str(float(self.rexi_ci_sy))
                                if self.rexi_ci_mu != None:
                                    idstr += '_mu'+str(float(self.rexi_ci_mu))
                            idstr += '_pr'+str(self.rexi_ci_primitive)


                    #idstr += '_rexithreadpar'+str(1 if self.rexi_thread_par else 0)

        return idstr


    def getRuntimeOptions(self):
        retval = ''

        if self.rexi_method != '' and self.rexi_method != None:
            retval += ' --rexi-method='+str(self.rexi_method)

            if self.rexi_method == 'direct':
                retval += ' --exp-direct-precompute-phin='+str(self.exp_direct_precompute_phin)

            else:
                retval += ' --rexi-sphere-preallocation='+str(self.rexi_sphere_preallocation)

                if self.rexi_method == 'file':

                    if self.p_job_dirpath == None:
                        raise Exception("self.p_job_dirpath not set!")

                    # REXI Files
                    file_params = []
                    for coeffs in self.rexi_files_coefficients:
                        coeff_filepath = self.p_job_dirpath+'/'+coeffs.unique_id_string
                        coeffs.write_file(coeff_filepath)
                        file_params.append(coeffs.function_name+':'+coeff_filepath)

                    retval += ' --rexi-files='+(",".join(file_params))

                elif self.rexi_method == 'terry':

                    # REXI Terry
                    retval += ' --rexi-terry-m='+str(self.rexi_terry_m)
                    retval += ' --rexi-terry-h='+str(self.rexi_terry_h)
                    retval += ' --rexi-terry-reduce-to-half='+str(self.rexi_terry_reduce_to_half)

                elif self.rexi_method == 'ci':

                    retval += ' --rexi-ci-primitive='+str(self.rexi_ci_primitive)
                    retval += ' --rexi-ci-n='+str(self.rexi_ci_n)

                    if self.rexi_ci_max_real != None:
                        retval += ' --rexi-ci-max-real='+str(self.rexi_ci_max_real)
                        if self.rexi_ci_max_imag != None:
                            retval += ' --rexi-ci-max-imag='+str(self.rexi_ci_max_imag)
                    else:
                        if self.rexi_ci_sx != None:
                            retval += ' --rexi-ci-sx='+str(self.rexi_ci_sx)
                        if self.rexi_ci_sy != None:
                            retval += ' --rexi-ci-sy='+str(self.rexi_ci_sy)
                        if self.rexi_ci_mu != None:
                            retval += ' --rexi-ci-mu='+str(self.rexi_ci_mu)
        return retval
    
    
    