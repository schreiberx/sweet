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

        self.gui = None

        self.polvani_rossby = None
        self.polvani_froude = None


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

        self.reuse_plans = -1
        self.comma_separated_tags = None

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

            if self.plane_domain_size == None:
                idstr += '_spap'+str(self.semi_lagrangian_approximate_sphere_geometry)


            if self.space_use_spectral_basis_diffs != 1:
                idstr += '_spd'+str(self.space_use_spectral_basis_diffs)

        if not 'runtime.reuse_plans' in filter_list:
            if self.reuse_plans != -1:
                idstr += '_plans'+str(self.reuse_plans)

        if not 'runtime.comma_separated_tags' in filter_list:
            if self.comma_separated_tags != None:
                idstr += '_tags'+str(self.comma_separated_tags)

        if idstr != '':
            idstr = "RT"+idstr

        for key, param in self.user_defined_parameters.items():
            idstr += '_'+param['id']+str(param['value'])

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

            if self.rexi_method != 'direct':
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

        retval += ' --semi-lagrangian-approximate-sphere-geometry='+str(self.semi_lagrangian_approximate_sphere_geometry)

        retval += ' --compute-error='+str(self.compute_error)

        retval += ' --reuse-plans='+str(self.reuse_plans)

        if self.comma_separated_tags != None:
            retval += ' --comma-separated-tags='+str(self.comma_separated_tags)


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

