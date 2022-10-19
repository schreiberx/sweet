#! /usr/bin/env python3

import sys

from mule_local.JobCompileOptions import *
from mule.InfoError import *

__all__ = ['JobRuntimeOptions']

class JobRuntimeOptions(InfoError):


    def __init__(self, dummy_init = False):

        self.init_phase = True

        InfoError.__init__(self, "JobRuntimeOptions")

        # String to job directory.
        # This is required to generate additional files such as REXI coefficients stored separately in each job directory
        self.p_job_dirpath = None


        self.space_res_spectral = None
        self.space_res_physical = None

        self.output_timestep_size = None
        self.output_filename = ''
        self.output_file_mode = ''

        self.f_sphere = None
        self.verbosity = 0

        self.instability_checks = None

        # Use 14 digits per default
        self.floating_point_output_digits = 12

        self.timestepping_method = None
        self.timestepping_order = 1
        self.timestepping_order2 = 1

        self.semi_lagrangian_max_iterations = None
        self.semi_lagrangian_convergence_threshold = None

        self.timestep_size = None
        self.max_timesteps_nr = -1

        self.normal_mode_analysis = None

        #
        # REXI method:
        #
        # 'file': default & can be used to also load terry, cauchy & butcher coefficients via text files
        # 'direct': Use direct solution, if available
        # 'terry': deprecated T-REXI method
        # 'butcher': deprecated Butcher-REXI method
        #
        self.rexi_method = 'file'

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

        ## Parameters for direct exp. solver
        self.exp_direct_precompute_phin = 0;

        self.gui = None

        self.polvani_rossby = None
        self.polvani_froude = None

        # Parameters for LibPFASST
        self.libpfasst_nlevels = None
        self.libpfasst_niters = None
        self.libpfasst_nnodes = None
        self.libpfasst_nsweeps = None
        self.libpfasst_nodes_type = None
        self.libpfasst_coarsening_multiplier = None
        self.libpfasst_use_rexi = None
        self.libpfasst_implicit_coriolis_force = None
        self.libpfasst_use_rk_stepper = None
        self.libpfasst_u2 = None
        self.libpfasst_u4 = None
        self.libpfasst_u6 = None
        self.libpfasst_u8 = None
        self.libpfasst_u_fields = None

        self.gravitation= None
        self.h0 = None
        self.sphere_rotating_coriolis_omega = None
        self.sphere_radius = None

        self.plane_domain_size = None    # Plane: Domain size

        # Specify benchmark name
        self.benchmark_name = None

        self.benchmark_galewsky_umax = -1
        self.benchmark_galewsky_hamp = -1
        self.benchmark_galewsky_phi2 = -1

        self.benchmark_normal_modes_case = None

        self.semi_lagrangian_approximate_sphere_geometry = 0

        self.space_grid_use_c_staggering = 0
        self.space_use_spectral_basis_diffs = 1
        self.viscosity = None
        self.viscosity_order = None

        #self.uselineardiv = None
        self.use_nonlinear_only_visc = None
        self.advection_rotation_angle = None
        self.advection_velocity = None
        self.max_simulation_time = 0.001
        self.max_wallclock_time = -1

        self.compute_error = 0

        self.reuse_plans = "quick"
        self.comma_separated_tags = None

        ## parareal parameters
        self.parareal_enabled = 0
        self.parareal_coarse_slices = None
        self.parareal_convergence_threshold = -1
        self.parareal_verbosity = 0
        self.parareal_max_simulation_time = None
        self.parareal_coarse_timestepping_method = None
        self.parareal_coarse_timestepping_order = 1
        self.parareal_coarse_timestepping_order2 = 1
        self.parareal_coarse_timestep_size = -1;
        self.parareal_load_ref_csv_files = 0;
        self.parareal_path_ref_csv_files = "";
        self.parareal_load_fine_csv_files = 0;
        self.parareal_path_fine_csv_files = "";
        self.parareal_store_iterations = 1;
        self.parareal_spatial_coarsening = None;
        self.parareal_max_iter = None;


        ## XBraid parameters
        self.xbraid_enabled = 0;
        self.xbraid_max_levels = None;
        self.xbraid_skip = None;
        self.xbraid_min_coarse = None;
        self.xbraid_nrelax = None;
        self.xbraid_nrelax0 = None;
        self.xbraid_tol = None;
        self.xbraid_tnorm = None;
        self.xbraid_cfactor = None;
        self.xbraid_cfactor0 = None;
        self.xbraid_max_iter = None;
        self.xbraid_fmg = None;
        self.xbraid_fmg_vcyc = None;
        self.xbraid_res = None;
        self.xbraid_storage = None;
        self.xbraid_print_level = None;
        self.xbraid_access_level = None;
        self.xbraid_run_wrapper_tests = None;
        self.xbraid_fullrnorm = None;
        self.xbraid_use_seq_soln = None;
        self.xbraid_use_rand = None;
        self.xbraid_pt = None;
        self.xbraid_timestepping_method = None;
        self.xbraid_timestepping_order = None;
        self.xbraid_timestepping_order2 = None;
        self.xbraid_verbosity = None;
        self.xbraid_load_ref_csv_files = None;
        self.xbraid_path_ref_csv_files = None;
        self.xbraid_load_fine_csv_files = None;
        self.xbraid_path_fine_csv_files = None;
        self.xbraid_store_iterations = None;
        self.xbraid_spatial_coarsening = None;


        ## ODE parameters
        self.function_param_y0_real = None;
        self.function_param_y0_imag = None;
        self.function_param_L = None;
        self.function_param_N = None;
        self.function_param_extra = None;
        self.ode_model = None;

        #
        # User defined parameters
        # Each new entry must set three values:
        # 
        # 'id': The id for the unique id
        # 'value': The value of the parameter
        # 'option': The program option
        #
        # example:
        # self.user_defined_parameters['test-mode'] = {'id': 'tm', 'value': 1, 'option': '--test-mode='}
        #
        self.user_defined_parameters = {}

        self.init_phase = False


    def __setattr__(self, name, value):

        if name != 'init_phase':
            if not self.init_phase:
                if not name in self.__dict__:
                    raise Exception("Attribute '"+name+"' does not exist!")

        self.__dict__[name] = value



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

        if 'polvani_rossby' in d:
            self.polvani_rossby = float(d['polvani_rossby'])

        if 'polvani_froude' in d:
            self.polvani_froude = float(d['polvani_froude'])

        if 'timestep_size' in d:
            self.timestep_size = float(d['timestep_size'])

        if 'reuse_plans' in d:
            self.reuse_plans = int(d['reuse_plans'])

        if 'comma_separated_tags' in d:
            self.comma_separated_tags = d['comma_separated_tags']

        if 'exp_direct_precompute_phin' in d:
            self.exp_direct_precompute_phin = d['exp_direct_precompute_phin']



    def getUniqueID(self, compileOptions : JobCompileOptions, filter_list : list = []):
        idstr = ''

        if not 'runtime.benchmark' in filter_list:
            if self.benchmark_name != None:
                idstr += '_b'+str(self.benchmark_name)

        if not 'runtime.galewsky_params' in filter_list:
            if self.benchmark_galewsky_umax > 0:
                idstr += '_bgu'+str("{:.4E}".format(self.benchmark_galewsky_umax))

            if self.benchmark_galewsky_hamp > 0:
                idstr += '_bgh'+str("{:.4E}".format(self.benchmark_galewsky_hamp))

            if self.benchmark_galewsky_phi2 > 0:
                idstr += '_bgp'+str("{:.4E}".format(self.benchmark_galewsky_phi2))

        if not 'runtime.normal_modes_params' in filter_list:
            if self.benchmark_normal_modes_case != None:
                idstr += '_bcase'+str(self.benchmark_normal_modes_case)

        if not 'runtime.simparams' in filter_list:
            if self.gravitation!= None:
                idstr += '_g'+str("{:05.2f}".format(self.gravitation))
            if self.h0 != None:
                idstr += '_h'+str("{:010.3f}".format(self.h0))
            if self.sphere_rotating_coriolis_omega != None:
                idstr += '_f'+str("{:e}".format(self.sphere_rotating_coriolis_omega))

            if compileOptions.sphere_spectral_space == 'enable':
                if self.sphere_radius != None:
                    idstr += '_a'+str(self.sphere_radius)
                if self.f_sphere != None:
                    idstr += '_fsph'+str(self.f_sphere)

            if self.viscosity != None:
                idstr += '_u'+str(self.viscosity)
            if self.viscosity_order != None:
                idstr += '_U'+str(self.viscosity_order)

            if self.advection_rotation_angle != None:
                idstr += '_ar'+str(self.advection_rotation_angle)

            if self.advection_velocity != None:
                idstr += '_av'+str(self.advection_velocity).replace(",", "_")

        if 'timestep' in filter_list:
            raise Exception("Deprecated")

        if not 'runtime.timestepping' in filter_list:
            if self.timestepping_method != None:
                if not 'runtime.timestepping_method' in filter_list:
                    idstr += '_tsm_'+self.timestepping_method

                if not 'runtime.timestepping_order' in filter_list:
                    idstr += '_tso'+str(self.timestepping_order)
                    idstr += '_tsob'+str(self.timestepping_order2)

                if not 'runtime.semi_lagrangian' in filter_list:
                    if self.semi_lagrangian_max_iterations != None:
                        idstr += '_sli'+str(self.semi_lagrangian_max_iterations)

                    if self.semi_lagrangian_convergence_threshold != None:
                        idstr += '_slc'+str("{:0.5e}".format(self.semi_lagrangian_convergence_threshold))


            if not 'runtime.timestepping_size' in filter_list:
                if self.timestep_size != None:
                    # Leading number is the total number of digits!
                    if self.timestep_size < 1:
                        idstr += '_dt'+str("{:0.8e}".format(self.timestep_size))
                    else:
                        idstr += '_dt'+str("{:08.2f}".format(self.timestep_size))

            if not 'runtime.max_timesteps_nr' in filter_list:
                if self.max_timesteps_nr != -1:
                    idstr += '_T'+str(self.max_timesteps_nr).zfill(3)

            if not 'runtime.max_wallclock_time' in filter_list:
                idstr += '_W'+str(self.max_wallclock_time).zfill(6)


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


        if not 'runtime.polvani' in filter_list:
            if self.polvani_rossby != None:
                idstr += '_PR'+str(self.polvani_rossby)

            if self.polvani_froude != None:
                idstr += '_PF'+str(self.polvani_froude)

        # LibPFASST attributes
        if not 'runtime.libpfasst' in filter_list:
            if self.libpfasst_nlevels != None:
                idstr += '_pf_nlev'+str(self.libpfasst_nlevels)
            if self.libpfasst_niters != None:
                idstr += '_pf_nit'+str(self.libpfasst_niters)
            if self.libpfasst_nnodes != None:
                idstr += '_pf_nnod'+str(self.libpfasst_nnodes)
            if self.libpfasst_nsweeps != None:
                idstr += '_pf_nswps'+str(self.libpfasst_nsweeps)
            if self.libpfasst_nodes_type != None:
                idstr += '_'+str(self.libpfasst_nodes_type)
            if self.libpfasst_coarsening_multiplier != None:
                idstr += '_pf_coars'+str(self.libpfasst_coarsening_multiplier)
            if self.libpfasst_use_rexi != None:
                idstr += '_pf_rexi'+str(self.libpfasst_use_rexi)
            if self.libpfasst_implicit_coriolis_force != None:
                idstr += '_pf_icf'+str(self.libpfasst_implicit_coriolis_force)
            if self.libpfasst_use_rk_stepper != None:
                idstr += '_pf_rk'+str(self.libpfasst_use_rk_stepper)

        if not 'runtime.disc_space' in filter_list:
            if self.space_res_spectral != None:
                if isinstance(self.space_res_spectral, (list, tuple)):
                    idstr += '_M'+str("x".join([str(x).zfill(4) for x in self.space_res_spectral]))
                else:
                    idstr += '_M'+str(self.space_res_spectral).zfill(4)

            if self.space_res_physical != None:
                if isinstance(self.space_res_physical, (list, tuple)):
                    idstr += '_N'+str("x".join([str(x).zfill(4) for x in self.space_res_physical]))
                else:
                    idstr += '_N'+str(self.space_res_physical).zfill(4)

            if self.plane_domain_size != None:
                if isinstance(self.plane_domain_size, (list, tuple)):
                    idstr += '_X'+str("x".join([str(x).zfill(4) for x in self.plane_domain_size]))
                else:
                    idstr += '_X'+str(self.plane_domain_size)

            if self.semi_lagrangian_approximate_sphere_geometry == None:
                idstr += '_spap'+str(self.semi_lagrangian_approximate_sphere_geometry)


            if self.space_use_spectral_basis_diffs != 1:
                idstr += '_spd'+str(self.space_use_spectral_basis_diffs)

        if not 'runtime.reuse_plans' in filter_list:
            if self.reuse_plans != -1:
                idstr += '_plans'+str(self.reuse_plans)

        if not 'runtime.comma_separated_tags' in filter_list:
            if self.comma_separated_tags != None:
                idstr += '_tags'+str(self.comma_separated_tags)

        if not 'runtime.parareal' in filter_list:
            if self.parareal_enabled:
                if not 'runtime.parareal_coarse_slices' in filter_list:
                    idstr += '_par_'+str(self.parareal_coarse_slices)
                if not 'runtime.parareal_coarse_timestepping_method' in filter_list:
                    idstr += '_ptsm_'+str(self.parareal_coarse_timestepping_method)
                if not 'runtime.parareal_coarse_timestep_size' in filter_list:
                    idstr += '_pDt_'+str(self.parareal_coarse_timestep_size)
                if not 'runtime.parareal_store_iterations' in filter_list:
                    idstr += '_pStore_'+str(self.parareal_store_iterations)
            if not 'runtime.parareal_spatial_coarsening' in filter_list:
                    idstr += '_pSpc_'+str(self.parareal_spatial_coarsening)
            if not 'runtime.parareal_max_iter' in filter_list:
                    idstr += '_pMaxIter_'+str(self.parareal_max_iter)


        if not 'runtime.xbraid' in filter_list:
            if self.xbraid_enabled:
                if self.xbraid_max_levels != None:
                    idstr += '_xb_max_l'+str(self.xbraid_max_levels)
                if self.xbraid_skip!= None:
                    idstr += '_xb_skip'+str(self.xbraid_skip)
                if self.xbraid_min_coarse != None:
                    idstr += '_xb_min_c'+str(self.xbraid_min_coarse)
                if self.xbraid_nrelax != None:
                    idstr += '_xb_nrlx'+str(self.xbraid_nrelax)
                if self.xbraid_nrelax0 != None:
                    idstr += '_xb_nrlx0'+str(self.xbraid_nrelax0)
                if self.xbraid_tol != None:
                    idstr += '_xb_tol'+str(self.xbraid_tol)
                if self.xbraid_tnorm != None:
                    idstr += '_xb_tnorm'+str(self.xbraid_tnorm)
                if self.xbraid_cfactor != None:
                    idstr += '_xb_cfr'+str(self.xbraid_cfactor)
                if self.xbraid_cfactor0 != None:
                    idstr += '_xb_cfr0'+str(self.xbraid_cfactor0)
                if self.xbraid_max_iter != None:
                    idstr += '_xb_max_i'+str(self.xbraid_max_iter)
                if self.xbraid_fmg != None:
                    idstr += '_xb_fmg'+str(self.xbraid_fmg)
                if self.xbraid_fmg_vcyc != None:
                    idstr += '_xb_fmg_vcyc'+str(self.xbraid_fmg_vcyc)
                if self.xbraid_res != None:
                    idstr += '_xb_res'+str(self.xbraid_res)
                if self.xbraid_storage != None:
                    idstr += '_xb_stg'+str(self.xbraid_storage)
                if self.xbraid_print_level != None:
                    idstr += '_xb_prt_l'+str(self.xbraid_print_level)
                if self.xbraid_access_level != None:
                    idstr += '_xb_acc_l'+str(self.xbraid_access_level)
                if self.xbraid_run_wrapper_tests != None:
                    idstr += '_xb_wrap'+str(self.xbraid_run_wrapper_tests)
                if self.xbraid_fullrnorm != None:
                    idstr += '_xb_frnm'+str(self.xbraid_fullrnorm)
                if self.xbraid_use_seq_soln != None:
                    idstr += '_xb_seq'+str(self.xbraid_use_seq_soln)
                if self.xbraid_use_rand != None:
                    idstr += '_xb_rand'+str(self.xbraid_use_rand)
                if self.xbraid_pt != None:
                    idstr += '_xb_pt'+str(self.xbraid_pt)
                if self.xbraid_timestepping_method != None:
                    idstr += '_xb_tsm'+str(self.xbraid_timestepping_method)
                if self.xbraid_timestepping_order != None:
                    idstr += '_xb_tso'+str(self.xbraid_timestepping_order)
                if self.xbraid_timestepping_order != None:
                    idstr += '_xb_tso2'+str(self.xbraid_timestepping_order2)
                if self.xbraid_verbosity != None:
                    idstr += '_xb_verb'+str(self.xbraid_verbosity)
                if self.xbraid_load_ref_csv_files != None:
                    idstr += '_xb_load_ref'+str(self.xbraid_load_ref_csv_files)
                if self.xbraid_path_ref_csv_files != None:
                    idstr += '_xb_path_ref'+str(self.xbraid_path_ref_csv_files)
                if self.xbraid_load_fine_csv_files != None:
                    idstr += '_xb_load_fine'+str(self.xbraid_load_fine_csv_files)
                if self.xbraid_path_fine_csv_files != None:
                    idstr += '_xb_path_fine'+str(self.xbraid_path_fine_csv_files)
                if self.xbraid_store_iterations != None:
                    idstr += '_xb_store_iterations'+str(self.xbraid_store_iterations)
                if self.xbraid_spatial_coarsening != None:
                    idstr += '_xb_spc'+str(self.xbraid_spatial_coarsening)

        if not 'runtime.ode' in filter_list:
            if self.function_param_y0_real != None:
                idstr += '_ode_y0Re'+str(self.function_param_y0_real)
            if self.function_param_y0_imag != None:
                idstr += '_ode_y0Im'+str(self.function_param_y0_imag)
            if self.function_param_L != None:
                idstr += '_ode_L'+str(self.function_param_L)
            if self.function_param_N != None:
                idstr += '_ode_N'+str(self.function_param_N)
            if self.function_param_extra != None:
                idstr += '_ode_extra'+str(self.function_param_extra)
            if self.ode_model != None:
                idstr += '_ode_model'+str(self.ode_model)

        if idstr != '':
            idstr = "RT"+idstr

        for key, param in self.user_defined_parameters.items():
            idstr += '_'+param['id']+str(param['value'])

        print(idstr)

        return idstr


    def getRuntimeOptions(self):
        retval = ''

        self.cleanup_options()

        if self.gui != None:
            retval += ' -G '+str(self.gui)

        if self.gravitation!= None:
            retval += ' -g '+str(self.gravitation)
        if self.h0 != None:
            retval += ' -H '+str(self.h0)
        if self.sphere_rotating_coriolis_omega != None:
            retval += ' -f '+str(self.sphere_rotating_coriolis_omega)
        if self.sphere_radius != None:
            retval += ' -a '+str(self.sphere_radius)
        if self.f_sphere != None:
            retval += ' -F '+str(self.f_sphere)

        if self.space_res_spectral != None:
            if isinstance(self.space_res_spectral, (list, tuple)):
                retval += ' -M '+str(",".join([str(x) for x in self.space_res_spectral]))
            else:
                retval += ' -M '+str(self.space_res_spectral)

        if self.space_res_physical != None:
            if isinstance(self.space_res_physical, (list, tuple)):
                retval += ' -N '+str(",".join([str(x) for x in self.space_res_physical]))
            else:
                retval += ' -N '+str(self.space_res_physical)

        retval += ' --space-grid-use-c-staggering='+str(self.space_grid_use_c_staggering)
        retval += ' -S '+str(self.space_use_spectral_basis_diffs)

        if self.plane_domain_size != None:
            if isinstance(self.plane_domain_size, (int, float)):
                retval += ' -X '+str(self.plane_domain_size)
                retval += ' -Y '+str(self.plane_domain_size)
            else:
                retval += ' -X '+str(self.plane_domain_size[0])
                retval += ' -Y '+str(self.plane_domain_size[1])

        if self.benchmark_name != None:
            retval += ' --benchmark-name='+str(self.benchmark_name)

        if self.benchmark_normal_modes_case != None:
            retval += ' --benchmark-normal-modes-case='+str(self.benchmark_normal_modes_case)


        retval += ' -v '+str(self.verbosity)

        if self.timestep_size != None:
            retval += ' --dt='+str(self.timestep_size)

        if self.max_timesteps_nr != -1:
            retval += ' -T '+str(self.max_timesteps_nr)

        if self.output_timestep_size != None:
            retval += ' -o '+str(self.output_timestep_size)

        if self.output_filename != '':
            retval += ' --output-file-name='+self.output_filename
        elif self.output_timestep_size == None:
            retval += ' --output-file-name=-'

        if self.output_file_mode != '':
            retval += ' --output-file-mode='+self.output_file_mode


        if self.viscosity != None:
            retval += ' -u '+str(self.viscosity)

        if self.viscosity_order != None:
            retval += ' -U '+str(self.viscosity_order)

        retval += ' -t '+str(self.max_simulation_time)

        retval += ' --max-wallclock-time '+str(self.max_wallclock_time)

        if self.instability_checks != None:
            retval += ' --instability-checks='+str(self.instability_checks)

        if self.floating_point_output_digits >= 0:
            retval += ' -d '+str(self.floating_point_output_digits)

        #if self.uselineardiv != None:
        #    retval += ' --use-only-linear-div='+str(self.uselineardiv)

        if self.use_nonlinear_only_visc != None:
            retval += ' --use-nonlinear-only-visc='+str(self.use_nonlinear_only_visc)

        if self.advection_rotation_angle != None:
            retval += ' --advection-rotation-angle='+str(self.advection_rotation_angle)

        if self.advection_velocity != None:
            retval += ' --advection-velocity='+str(self.advection_velocity)

        if self.timestepping_method != None:
            retval += ' --timestepping-method='+self.timestepping_method
            retval += ' --timestepping-order='+str(self.timestepping_order)
            retval += ' --timestepping-order2='+str(self.timestepping_order2)


        if self.semi_lagrangian_max_iterations != None:
            retval += ' --semi-lagrangian-max-iterations='+str(self.semi_lagrangian_max_iterations)

        if self.semi_lagrangian_convergence_threshold != None:
            retval += ' --semi-lagrangian-convergence-threshold='+str(self.semi_lagrangian_convergence_threshold)

        if self.normal_mode_analysis != None:
            retval += ' --normal-mode-analysis-generation='+str(self.normal_mode_analysis)

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


        if self.polvani_rossby != None:
            retval += ' --polvani-rossby='+str(self.polvani_rossby)

        if self.polvani_froude != None:
            retval += ' --polvani-froude='+str(self.polvani_froude)

        # TODO: think about a way to check that we're in LibPFASST mode s.th. we do not have to check every if statement
        if self.libpfasst_nlevels != None:
            retval += ' --libpfasst-nlevels='+str(self.libpfasst_nlevels)
        if self.libpfasst_niters != None:
            retval += ' --libpfasst-niters='+str(self.libpfasst_niters)
        if self.libpfasst_nnodes != None:
            retval += ' --libpfasst-nnodes='+str(self.libpfasst_nnodes)
        if self.libpfasst_nsweeps != None:
            retval += ' --libpfasst-nsweeps='+str(self.libpfasst_nsweeps)
        if self.libpfasst_nodes_type != None:
            retval += ' --libpfasst-nodes-type='+str(self.libpfasst_nodes_type)
        if self.libpfasst_coarsening_multiplier != None:
            retval += ' --libpfasst-coarsening-multiplier='+str(self.libpfasst_coarsening_multiplier)
        if self.libpfasst_use_rexi != None:
            retval += ' --libpfasst-use-rexi='+str(self.libpfasst_use_rexi)
        if self.libpfasst_implicit_coriolis_force != None:
            retval += ' --libpfasst-implicit-coriolis-force='+str(self.libpfasst_implicit_coriolis_force)
        if self.libpfasst_use_rk_stepper != None:
            retval += ' --libpfasst-use-rk-stepper='+str(self.libpfasst_use_rk_stepper)
        if self.libpfasst_u2 != None:
            retval += ' --libpfasst-u2='+str(self.libpfasst_u2)
        if self.libpfasst_u4 != None:
            retval += ' --libpfasst-u4='+str(self.libpfasst_u4)
        if self.libpfasst_u6 != None:
            retval += ' --libpfasst-u6='+str(self.libpfasst_u6)
        if self.libpfasst_u8 != None:
            retval += ' --libpfasst-u8='+str(self.libpfasst_u8)
        if self.libpfasst_u_fields != None:
            retval += ' --libpfasst-u-fields='+str(self.libpfasst_u_fields)

        retval += ' --semi-lagrangian-approximate-sphere-geometry='+str(self.semi_lagrangian_approximate_sphere_geometry)

        retval += ' --compute-error='+str(self.compute_error)

        retval += ' --reuse-plans='+str(self.reuse_plans)

        if self.comma_separated_tags != None:
            retval += ' --comma-separated-tags='+str(self.comma_separated_tags)

        ## Parareal parameters
        if self.parareal_enabled:
            retval += " --parareal-enable=1"
            retval += " --parareal-coarse-slices="+str(self.parareal_coarse_slices)
            retval += " --parareal-convergence-threshold="+str(self.parareal_convergence_threshold)
            retval += " --parareal-verbosity="+str(self.parareal_verbosity)
            retval += " --parareal-max-simulation-time="+str(self.parareal_max_simulation_time)
            retval += " --parareal-coarse-timestepping-method="+str(self.parareal_coarse_timestepping_method)
            retval += " --parareal-coarse-timestepping-order="+str(self.parareal_coarse_timestepping_order)
            retval += " --parareal-coarse-timestepping-order2="+str(self.parareal_coarse_timestepping_order2)
            retval += " --parareal-coarse-timestep-size="+str(self.parareal_coarse_timestep_size);
            retval += " --parareal-load-ref-csv-files="+str(self.parareal_load_ref_csv_files);
            retval += " --parareal-path-ref-csv-files="+str(self.parareal_path_ref_csv_files);
            retval += " --parareal-load-fine-csv-files="+str(self.parareal_load_fine_csv_files);
            retval += " --parareal-path-fine-csv-files="+str(self.parareal_path_fine_csv_files);
            retval += " --parareal-store-iterations="+str(self.parareal_store_iterations);
            retval += " --parareal-spatial-coarsening="+str(self.parareal_spatial_coarsening);
            if self.parareal_max_iter != None:
                retval += " --parareal-max-iter="+str(self.parareal_max_iter);

            ##if self.parareal_coarse_timestep_size > 0:
            ##    retval += " --parareal-coarse-timestep-size="+str(self.parareal_coarse_timestep_size);
            ##else:
            ##    retval += " --parareal-coarse-timestep-size="+str(self.parareal_max_simulation_time/self.parareal_coarse_slices);

        ## XBraid parameters
        if self.xbraid_enabled:
            retval += " --xbraid-enable=1"
            retval += " --xbraid-max-levels="+str(self.xbraid_max_levels)
            retval += " --xbraid-skip="+str(self.xbraid_skip)
            retval += " --xbraid-min-coarse="+str(self.xbraid_min_coarse)
            retval += " --xbraid-nrelax="+str(self.xbraid_nrelax)
            retval += " --xbraid-nrelax0="+str(self.xbraid_nrelax0)
            retval += " --xbraid-tol="+str(self.xbraid_tol)
            retval += " --xbraid-tnorm="+str(self.xbraid_tnorm)
            retval += " --xbraid-cfactor="+str(self.xbraid_cfactor)
            retval += " --xbraid-cfactor0="+str(self.xbraid_cfactor0)
            retval += " --xbraid-max-iter="+str(self.xbraid_max_iter)
            retval += " --xbraid-fmg="+str(self.xbraid_fmg)
            retval += " --xbraid-fmg-vcyc="+str(self.xbraid_fmg_vcyc)
            retval += " --xbraid-res="+str(self.xbraid_res)
            retval += " --xbraid-storage="+str(self.xbraid_storage)
            retval += " --xbraid-print-level="+str(self.xbraid_print_level)
            retval += " --xbraid-access-level="+str(self.xbraid_access_level)
            retval += " --xbraid-run-wrapper-tests="+str(self.xbraid_run_wrapper_tests)
            retval += " --xbraid-fullrnorm="+str(self.xbraid_fullrnorm)
            retval += " --xbraid-use-seq-soln="+str(self.xbraid_use_seq_soln)
            retval += " --xbraid-use-rand="+str(self.xbraid_use_rand)
            retval += " --xbraid-pt="+str(self.xbraid_pt)
            retval += " --xbraid-timestepping-method="+str(self.xbraid_timestepping_method)
            retval += " --xbraid-timestepping-order="+str(self.xbraid_timestepping_order2)
            retval += " --xbraid-timestepping-order2="+str(self.xbraid_timestepping_order2)
            retval += " --xbraid-verbosity="+str(self.xbraid_verbosity)
            retval += " --xbraid-load-ref-csv-files="+str(self.xbraid_load_ref_csv_files)
            retval += " --xbraid-path-ref-csv-files="+str(self.xbraid_path_ref_csv_files)
            retval += " --xbraid-load-fine-csv-files="+str(self.xbraid_load_fine_csv_files)
            retval += " --xbraid-path-fine-csv-files="+str(self.xbraid_path_fine_csv_files)
            retval += " --xbraid-store-iterations="+str(self.xbraid_store_iterations)
            retval += " --xbraid-spatial-coarsening="+str(self.xbraid_spatial_coarsening)

        ## ODE parameters
        if self.function_param_y0_real != None:
            retval += ' --function-param-y0-real='+str(self.function_param_y0_real)
        if self.function_param_y0_imag != None:
            retval += ' --function-param-y0-imag='+str(self.function_param_y0_imag)
        if self.function_param_L != None:
            retval += ' --function-param-L='+str(self.function_param_L)
        if self.function_param_N != None:
            retval += ' --function-param-N='+str(self.function_param_N)
        if self.function_param_extra != None:
            retval += ' --function-param-extra='+str(self.function_param_extra)
        if self.ode_model != None:
            retval += ' --ode-model='+str(self.ode_model)


        for key, param in self.user_defined_parameters.items():
            retval += ' '+param['option']+str(param['value'])

        return retval


    def cleanup_options(self):

        # Cleanup old format of reuse plans to new one
        if self.reuse_plans == -1:
            self.reuse_plans = "quick"
        elif self.reuse_plans == 0:
            self.reuse_plans = "save"
        elif self.reuse_plans == 1:
            self.reuse_plans = "load"
        elif self.reuse_plans == 2:
            self.reuse_plans = "require_load"


    def get_jobscript_plan_exec_prefix(self, compile, runtime):
        content = ""

        runtime.cleanup_options()

        plan_files = []
        if compile.plane_spectral_space == 'enable':
            plan_files.append('sweet_fftw')

        if compile.sphere_spectral_space == 'enable':
            plan_files.append('shtns_cfg')
            plan_files.append('shtns_cfg_fftw')


        if runtime.reuse_plans == "quick":
            # Quick plan generation mode, nothing to do
            pass

        elif "load" in runtime.reuse_plans:
            content += "\n"
            # Reuse plans if available
            # No error if plans don't exist
            for i in plan_files:
                if runtime.reuse_plans == "require_load":
                    content += "cp ../"+i+" ./ || exit 1\n"
                else:
                    content += "cp ../"+i+" ./ 2>/dev/null\n"

        elif runtime.reuse_plans == "save":
            # Create plans, don't load/store them
            pass

        else:
            raise Exception("Invalid reuse_plans value '"+str(runtime.reuse_plans)+"'")

        return content



    def get_jobscript_plan_exec_suffix(self, compile, runtime):
        content = ""

        runtime.cleanup_options()

        plan_files = []
        if compile.plane_spectral_space == 'enable':
            plan_files.append('sweet_fftw')

        if compile.sphere_spectral_space == 'enable':
            plan_files.append('shtns_cfg')
            plan_files.append('shtns_cfg_fftw')

        #
        # Reusing plans assumes them to be stored in the folder one level up in the hierarchy
        #
        if runtime.reuse_plans == "quick":
            # Quick plan generation mode, nothing to do
            pass

        elif "save" in runtime.reuse_plans:
            # Write back plans to main directory to reuse them
            content += "\n"
            for i in plan_files:
                content += "cp ./"+i+" ../ 2>/dev/null\n"
            content += "\n"
            pass

        elif "load" in runtime.reuse_plans:
            pass

        else:
            raise Exception("Invalid reuse_plans value '"+str(runtime.reuse_plans)+"'")


        return content


if __name__ == "__main__":

    p = JobRuntimeOptions()
    i = p.getRuntimeOptions()
    p.info(i)

    p.info("FIN")

