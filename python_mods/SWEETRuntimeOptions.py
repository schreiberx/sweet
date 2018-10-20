#! /usr/bin/env python3

from SWEETCompileOptions import *
from InfoError import *

__all__ = ['SWEETRuntimeOptions']

class SWEETRuntimeOptions(InfoError):

	"""
	class RuntimeOption:
		def __init__(self):
			self.variable = None
			self.name = 'rexi_use_direct_solution'

		def load_from_dict(self, d):
			if self.name in d:
				self.variable = d[self.name]

		def getUniqueID(self, filter_list = []):
			if self.variable:
				return '_REXIDIR'

		def getRuntimeOptions(self):
			return ' --rexi-use-direct-solution='+str(self.variable)
	"""



	def __init__(self, dummy_init = False):

		self.init_phase = True

		InfoError.__init__(self, "SWEETRuntimeOptions")

		self.mode_res = None
		self.phys_res = None

		self.output_timestep_size = 0.0001
		self.output_filename = ''
		self.output_file_mode = ''

		self.f_sphere = 0
		self.verbosity = 0

		# Deactivated per default for more performance
		self.instability_checks = 0

		# Use 14 digits per default
		self.floating_point_output_digits = 12

		self.timestepping_method = None#'ln_erk'
		self.timestepping_order = 1
		self.timestepping_order2 = 1

		self.timestep_size = 0.001
		self.max_timesteps = -1

		self.normal_mode_analysis = None

		self.rexi_method = ''

		self.rexi_beta_cutoff = 0

		self.rexi_file_n = 0
		self.rexi_file_h = 0
		self.rexi_file_test_abs = 0
		self.rexi_file_max_error = 0
		self.rexi_file_faf_dir = None

		self.rexi_ci_n = 128
		self.rexi_ci_max_real = None
		self.rexi_ci_max_imag = None
		self.rexi_ci_sx = None
		self.rexi_ci_sy = None
		self.rexi_ci_mu = None
		self.rexi_ci_primitive = 'circle'
		self.rexi_ci_gaussian_filter_scale = None
		self.rexi_ci_gaussian_filter_dt_norm = None
		self.rexi_ci_gaussian_filter_exp_N = None

		self.rexi_terry_m = 0
		self.rexi_terry_l = 11
		self.rexi_terry_h = 0.15

		self.rexi_half_poles = 0
		self.rexi_extended_modes = 0
		self.rexi_normalization = 0
		self.rexi_sphere_preallocation = 0
		self.rexi_use_direct_solution = 0

		self.polvani_rossby = -1.0
		self.polvani_froude = -1.0




		#self.g = 1		# gravity
		#self.h = 100000	# avg. height

		self.g = None
		self.h = None
		self.f = None
		self.r = None

#		self.g = 9.81		# gravity
#		self.h = 10000	# avg. height
#		self.f = 0.000072921	# \Omega coriolis effect
#		self.r = 6371220	# Sphere: Radius
		self.domain_size = None	# Plane: Domain size

		# 3: gaussian breaking dam
		# 4: geostrophic balance test case

		# Create error if bench_id is not specified!
		# This variable is DEPRECATED!!! Use benchmark_name instead!!!
		self.bench_id = None

		# Specify benchmark name
		self.benchmark_name = None

		self.benchmark_galewsky_umax = -1
		self.benchmark_galewsky_hamp = -1
		self.benchmark_galewsky_phi2 = -1

		self.use_robert_functions = 1

		self.pde_id = 0
		self.staggering = 0
		self.spectralderiv = 1
		self.viscosity = None
		self.viscosity_order = None
		self.uselineardiv = None
		self.uselocalvisc = None
		self.advection_velocity = None
		self.simtime = 0.001

		self.compute_error = 0

		self.reuse_plans = -1

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

		if 'half_poles' in d:
			self.rexi_half_poles = d['half_poles']

		if 'extended_modes' in d:
			self.rexi_extended_modes = d['extended_modes']

		if 'normalization' in d:
			self.rexi_normalization = d['normalization']

		if 'sphere_preallocation' in d:
			self.rexi_sphere_preallocation = d['sphere_preallocation']

		if 'use_direct_solution' in d:
			self.rexi_use_direct_solution = d['use_direct_solution']

		if 'terry_m' in d:
			self.rexi_terry_m = d['terry_m']

		if 'terry_h' in d:
			self.rexi_terry_h = d['terry_h']

		if 'file_n' in d:
			self.rexi_file_n = d['file_n']

		if 'file_h' in d:
			self.rexi_file_h = d['file_h']

		if 'file_test_abs' in d:
			self.rexi_file_test_abs = d['file_test_abs']

		if 'file_max_error' in d:
			self.rexi_file_max_error = d['file_max_error']

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

		if 'ci_gaussian_filter_scale' in d:
			self.rexi_ci_gaussian_filter_scale = float(d['ci_gaussian_filter_scale'])

		if 'ci_gaussian_filter_dt_norm' in d:
			self.rexi_ci_gaussian_filter_dt_norm = float(d['ci_gaussian_filter_dt_norm'])

		if 'ci_gaussian_filter_exp_N' in d:
			self.rexi_ci_gaussian_filter_exp_N = float(d['ci_gaussian_filter_exp_N'])

		if 'polvani_rossby' in d:
			self.polvani_rossby = float(d['polvani_rossby'])

		if 'polvani_froude' in d:
			self.polvani_froude = float(d['polvani_froude'])

		if 'timestep_size' in d:
			self.timestep_size = float(d['timestep_size'])

		if 'reuse_plans' in d:
			self.reuse_plans = int(d['reuse_plans'])




	def getUniqueID(self, compileOptions : SWEETCompileOptions, filter_list : list = []):
		idstr = ''

		if not 'runtime.benchmark' in filter_list:
			if self.benchmark_name != None:
				idstr += '_b'+str(self.benchmark_name)
			elif self.bench_id != None:
				idstr += '_b'+str(self.bench_id)

		if not 'runtime.galewsky_params' in filter_list:
			if self.benchmark_galewsky_umax > 0:
				idstr += '_bgu'+str("{:.4E}".format(self.benchmark_galewsky_umax))

			if self.benchmark_galewsky_hamp > 0:
				idstr += '_bgh'+str("{:.4E}".format(self.benchmark_galewsky_hamp))

			if self.benchmark_galewsky_phi2 > 0:
				idstr += '_bgp'+str("{:.4E}".format(self.benchmark_galewsky_phi2))

		if not 'runtime.simparams' in filter_list:
			if self.g != None:
				idstr += '_g'+str("{:05.2f}".format(self.g))
			if self.h != None:
				idstr += '_h'+str("{:010.3f}".format(self.h))
			if self.f != None:
				idstr += '_f'+str("{:e}".format(self.f))

			#idstr += '_p'+str(self.pde_id)

			if compileOptions.sphere_spectral_space == 'enable':
				if self.r != None:
					idstr += '_a'+str(self.r)
				idstr += '_fsph'+str(self.f_sphere)

			if self.viscosity != None:
				idstr += '_u'+str(self.viscosity)
			if self.viscosity_order != None:
				idstr += '_U'+str(self.viscosity_order)

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

			if not 'runtime.timestepping_size' in filter_list:
				idstr += '_C'+str("{:08.3f}".format(self.timestep_size))

			if self.max_timesteps != -1:
				idstr += '_T'+str(self.max_timesteps).zfill(3)

			idstr += '_S'+str(self.simtime).zfill(6)


		if not 'runtime.rexi' in filter_list:
			if self.rexi_method != '':
				if self.rexi_method == 'direct' or self.rexi_use_direct_solution:
					idstr += '_REXIDIR'
				else:
					if self.rexi_method == "file":
						idstr += '_REXIFIL'
						if not 'runtime.rexi_params' in filter_list:
							idstr += '_n'+str(self.rexi_file_n).zfill(8)
							idstr += '_h'+str(self.rexi_file_h)
							idstr += '_teabs'+str(self.rexi_file_test_abs).zfill(3)
							idstr += '_mxer'+str(self.rexi_file_max_error)

					elif self.rexi_method == "terry":
						idstr += '_REXITER'
						if not 'runtime.rexi_params' in filter_list:
							idstr += '_m'+str(self.rexi_terry_m).zfill(8)
							idstr += '_h'+str(self.rexi_terry_h)

					elif self.rexi_method == "ci":
						idstr += '_REXICI'

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

							if self.rexi_ci_gaussian_filter_scale != None:
								idstr += '_gfs'+str( "{:.4E}".format(self.rexi_ci_gaussian_filter_scale))
							if self.rexi_ci_gaussian_filter_dt_norm != None:
								idstr += '_gfd'+str( "{:.4E}".format(self.rexi_ci_gaussian_filter_dt_norm))
							if self.rexi_ci_gaussian_filter_exp_N != None:
								idstr += '_gfe'+str( "{:.4E}".format(self.rexi_ci_gaussian_filter_exp_N))

					if not 'runtime.rexi_params' in filter_list:
						idstr += '_nrm'+str(self.rexi_normalization)
						idstr += '_hlf'+str(self.rexi_half_poles)
						idstr += '_bf'+str(self.rexi_beta_cutoff)

					#if self.plane_or_sphere == 'sphere':
					#idstr += '_pre'+str(self.rexi_sphere_preallocation)
					idstr += '_ext'+str(self.rexi_extended_modes).zfill(2)

					#idstr += '_rexithreadpar'+str(1 if self.rexi_thread_par else 0)


		if not 'runtime.polvani' in filter_list:
			if self.polvani_rossby >= 0:
				idstr += '_PR'+str(self.polvani_rossby)

			if self.polvani_froude >= 0:
				idstr += '_PF'+str(self.polvani_froude)


		if not 'disc_space' in filter_list:
			if self.mode_res != None:
				if isinstance(self.mode_res, (list, tuple)):
					idstr += '_M'+str("x".join([str(x).zfill(4) for x in self.mode_res]))
				else:
					idstr += '_M'+str(self.mode_res).zfill(4)

			if self.phys_res != None:
				if isinstance(self.phys_res, (list, tuple)):
					idstr += '_N'+str("x".join([str(x).zfill(4) for x in self.phys_res]))
				else:
					idstr += '_N'+str(self.phys_res).zfill(4)

			idstr += '_rob'+str(self.use_robert_functions)

			if self.spectralderiv != 1:
				idstr += '_spd'+str(self.spectralderiv)

		if self.reuse_plans != -1:
			idstr += '_plans'+str(self.reuse_plans)

		if idstr != '':
			idstr = "RT"+idstr

		for key, param in self.user_defined_parameters.items():
			idstr += '_'+param['id']+str(param['value'])

		return idstr


	def getRuntimeOptions(self):
		retval = ''
		if self.g != None:
			retval += ' -g '+str(self.g)
		if self.h != None:
			retval += ' -H '+str(self.h)
		if self.f != None:
			retval += ' -f '+str(self.f)
		if self.r != None:
			retval += ' -a '+str(self.r)
		retval += ' -F '+str(self.f_sphere)
		if self.mode_res != None:
			if isinstance(self.mode_res, (list, tuple)):
				retval += ' -M '+str(",".join([str(x) for x in self.mode_res]))
			else:
				retval += ' -M '+str(self.mode_res)

		if self.phys_res != None:
			if isinstance(self.phys_res, (list, tuple)):
				retval += ' -N '+str(",".join([str(x) for x in self.phys_res]))
			else:
				retval += ' -N '+str(self.phys_res)

		retval += ' --pde-id '+str(self.pde_id)
		retval += ' --staggering='+str(self.staggering)
		retval += ' -S '+str(self.spectralderiv)

		if self.domain_size != None:
			retval += ' -X '+str(self.domain_size[0])
			retval += ' -Y '+str(self.domain_size[1])

		if self.bench_id != None:
			retval += ' -s '+str(self.bench_id)

		if self.benchmark_name != None:
			retval += ' --benchmark-name='+str(self.benchmark_name)

		retval += ' -v '+str(self.verbosity)

		retval += ' --dt='+str(self.timestep_size)

		if self.max_timesteps != -1:
			retval += ' -T '+str(self.max_timesteps)

		retval += ' -o '+str(self.output_timestep_size)

		if self.output_filename != '':
			retval += ' --output-file-name='+self.output_filename

		if self.output_file_mode != '':
			retval += ' --output-file-mode='+self.output_file_mode

		if self.output_timestep_size < 0:
			retval += ' --output-file-name=-'

		if self.viscosity != None:
			retval += ' -u '+str(self.viscosity)

		if self.viscosity_order != None:
			retval += ' -U '+str(self.viscosity_order)

		retval += ' -t '+str(self.simtime)

		retval += ' --instability-checks='+str(self.instability_checks)

		if self.floating_point_output_digits >= 0:
			retval += ' -d '+str(self.floating_point_output_digits)

		if self.uselineardiv != None:
			retval += ' --use-linear-div='+str(self.uselineardiv)

		if self.uselocalvisc != None:
			retval += ' --use-local-visc='+str(self.uselocalvisc)

		if self.advection_velocity != None:
			retval += ' --advection-velocity='+str(self.advection_velocity)

		if self.timestepping_method != None:
			retval += ' --timestepping-method='+self.timestepping_method
			retval += ' --timestepping-order='+str(self.timestepping_order)
			retval += ' --timestepping-order2='+str(self.timestepping_order2)

		if self.normal_mode_analysis != None:
			retval += ' --normal-mode-analysis-generation='+str(self.normal_mode_analysis)

		retval += ' --rexi-method='+str(self.rexi_method)

		if self.rexi_method != '':
			if self.rexi_method == 'direct' or self.rexi_use_direct_solution:
				if self.rexi_method == 'direct':
					self.rexi_use_direct_solution = 1
					retval += ' --rexi-use-direct-solution='+str(self.rexi_use_direct_solution)
			else:
				retval += ' --rexi-half='+str(self.rexi_half_poles)
				retval += ' --rexi-normalization='+str(self.rexi_normalization)
				retval += ' --rexi-sphere-preallocation='+str(self.rexi_sphere_preallocation)
				retval += ' --rexi-use-direct-solution='+str(self.rexi_use_direct_solution)
				retval += ' --rexi-ext-modes='+str(self.rexi_extended_modes)

				if self.rexi_beta_cutoff != 0:
					retval += ' --rexi-beta-cutoff='+str(self.rexi_beta_cutoff)

				if self.rexi_method == 'terry':
					# REXI Terry
					retval += ' --rexi-terry-m='+str(self.rexi_terry_m)
					retval += ' --rexi-terry-h='+str(self.rexi_terry_h)

				elif self.rexi_method == 'file':
					# REXI File
					retval += ' --rexi-file-n='+str(self.rexi_file_n)
					retval += ' --rexi-file-h='+str(self.rexi_file_h)
					retval += ' --rexi-file-test-abs='+str(self.rexi_file_test_abs)
					retval += ' --rexi-file-max-error='+str(self.rexi_file_max_error)
					if self.rexi_file_faf_dir != None:
						retval += ' --rexi-file-faf-dir='+str(self.rexi_file_faf_dir)

				elif self.rexi_method == 'ci':
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

					retval += ' --rexi-ci-primitive='+str(self.rexi_ci_primitive)
					retval += ' --rexi-ci-gaussian-filter-scale='+str(self.rexi_ci_gaussian_filter_scale)
					retval += ' --rexi-ci-gaussian-filter-dt-norm='+str(self.rexi_ci_gaussian_filter_dt_norm)
					retval += ' --rexi-ci-gaussian-filter-exp-N='+str(self.rexi_ci_gaussian_filter_exp_N)


		retval += ' --polvani-rossby='+str(self.polvani_rossby)
		retval += ' --polvani-froude='+str(self.polvani_froude)

		retval += ' --use-robert-functions='+str(self.use_robert_functions)

		retval += ' --compute-error='+str(self.compute_error)

		retval += ' --reuse-plans='+str(self.reuse_plans)

		for key, param in self.user_defined_parameters.items():
			retval += ' '+param['option']+str(param['value'])

		return retval


	def get_jobscript_plan_exec_prefix(self, compile, runtime):
		content = ""

		plan_files = []
		if compile.plane_spectral_space == 'enable':
			plan_files.append('sweet_fftw')

		if compile.sphere_spectral_space == 'enable':
			plan_files.append('shtns_cfg')
			plan_files.append('shtns_cfg_fftw')

		#
		# Reusing plans assumes them to be stored in the folder one level up in the hierarchy
		#
		if runtime.reuse_plans == -1:
			# Quick plan generation mode, nothing to do
			pass

		elif runtime.reuse_plans == 0:
			# Create plans, don't load/store them
			pass

		elif runtime.reuse_plans == 1:
			content += "\n"
			# Reuse plans if available
			# No error if plans don't exist
			for i in plan_files:
				content += "cp ../"+i+" ./ 2>/dev/null\n"
				
		elif runtime.reuse_plans == 2:
			content += "\n"
			# Reuse and trigger error if they are not available
			for i in plan_files:
				content += "cp ../"+i+" ./ || exit 1\n"
		else:
			raise Exception("Invalid reuse_plans value"+str(jg.runtime.reuse_plans))

		return content



	def get_jobscript_plan_exec_suffix(self, compile, runtime):
		content = ""

		plan_files = []
		if compile.plane_spectral_space == 'enable':
			plan_files.append('sweet_fftw')

		if compile.sphere_spectral_space == 'enable':
			plan_files.append('shtns_cfg')
			plan_files.append('shtns_cfg_fftw')

		#
		# Reusing plans assumes them to be stored in the folder one level up in the hierarchy
		#
		if runtime.reuse_plans == -1:
			# Quick plan generation mode, nothing to do
			pass

		elif runtime.reuse_plans == 0:
			# Create plans, don't load/store them
			pass

		elif runtime.reuse_plans == 1:
			content += "\n"
			# Write back plans to main directory to reuse them
			for i in plan_files:
				content += "cp ./"+i+" ../ 2>/dev/null\n"
				
		elif runtime.reuse_plans == 2:
			pass

		else:
			raise Exception("Invalid reuse_plans value"+str(jg.runtime.reuse_plans))


		return content


if __name__ == "__main__":

	p = SWEETRuntimeOptions()
	i = p.getRuntimeOptions()
	p.info(i)

	p.info("FIN")

