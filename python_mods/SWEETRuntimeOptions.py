#! /usr/bin/env python3


class SWEETRuntimeOptions():

	def __init__(self):
	#	self.mode_res = 32
		self.mode_res = -1
		self.phys_res = -1

		self.output_timestep_size = 0.0001
		self.output_filename = ''

		self.f_sphere = 0
		self.verbosity = 0

		self.stability_checks = 0

		self.floating_point_output_digits = -1

		self.timestepping_method = 'ln_erk'
		self.timestepping_order = 1
		self.timestepping_order2 = 1

		self.timestep_size = 0.001
		self.max_timesteps = -1

		self.normal_mode_analysis = 0


		self.rexi_method = ''

		self.rexi_beta_cutoff = 0

		self.rexi_file_n = 0
		self.rexi_file_h = 0
		self.rexi_file_test_abs = 0
		self.rexi_file_max_error = 0
		self.rexi_file_faf_dir = None

		self.rexi_ci_n = 128
		self.rexi_ci_max_real = -999
		self.rexi_ci_max_imag = -999
		self.rexi_ci_sx = 50
		self.rexi_ci_sy = 50
		self.rexi_ci_mu = 0
		self.rexi_ci_primitive = 'circle'
		self.rexi_ci_gaussian_filter_scale = 0
		self.rexi_ci_gaussian_filter_dt_norm = 0
		self.rexi_ci_gaussian_filter_exp_N = 0

		self.rexi_m = 0
		self.rexi_l = 11
		self.rexi_h = 0.15

		self.rexi_half_poles = 1
		self.rexi_extended_modes = 0
		self.rexi_normalization = 0
		self.rexi_sphere_preallocation = 0
		self.rexi_use_direct_solution = 0

		self.polvani_rossby = -1.0
		self.polvani_froude = -1.0




		#self.g = 1		# gravity
		#self.h = 100000	# avg. height

		self.g = 9.81		# gravity
		self.h = 10000	# avg. height

		self.f = 0.000072921	# \Omega coriolis effect
		#self.f = 2.0*f	# Coriolis effect for constant f (not multiplied with 2.0)

		self.r = 6371220	# Sphere: Radius
		self.domain_size = 1	# Plane: Domain size

		# 3: gaussian breaking dam
		# 4: geostrophic balance test case
		self.bench_id = 4
		self.benchmark_name = ""

		self.benchmark_galewsky_umax = -1
		self.benchmark_galewsky_hamp = -1
		self.benchmark_galewsky_phi2 = -1

		self.use_robert_functions = 1

		self.pde_id = 0
		self.staggering = 0
		self.spectralderiv = 1
		self.uselineardiv = 0
		self.viscosity = 0
		self.viscosity_order = 0

		self.simtime = 0.001

		self.compute_error = 0

		return



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
			self.rexi_m = d['terry_m']

		if 'terry_h' in d:
			self.rexi_h = d['terry_h']

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

	def getUniqueID(self, compileOptions):
		idstr = ''

		if False:
			if self.benchmark_name != '':
				idstr += '_b'+str(self.benchmark_name)
			else:
				idstr += '_b'+str(self.bench_id)

		if self.benchmark_galewsky_umax > 0:
			idstr += '_bgu'+str("{:.4E}".format(self.benchmark_galewsky_umax))

		if self.benchmark_galewsky_hamp > 0:
			idstr += '_bgh'+str("{:.4E}".format(self.benchmark_galewsky_hamp))

		if self.benchmark_galewsky_phi2 > 0:
			idstr += '_bgp'+str("{:.4E}".format(self.benchmark_galewsky_phi2))

		idstr += '_g'+str(self.g)
		idstr += '_h'+str(self.h)
		idstr += '_f'+str(self.f)

		#idstr += '_p'+str(self.pde_id)

		if compileOptions.plane_or_sphere == 'sphere':
			idstr += '_a'+str(self.r)
			idstr += '_u'+str(self.viscosity)
#			idstr += '_U'+str(self.viscosity_order)

			#idstr += '_rob'+str(self.use_robert_functions)
			idstr += '_fsph'+str(self.f_sphere)

#		idstr += '_t'+str(self.simtime).zfill(8)
#		idstr += '_o'+str(self.output_timestep_size).zfill(8)

		idstr += '_tsm_'+self.timestepping_method
		idstr += '_tso'+str(self.timestepping_order)
		idstr += '_tsob'+str(self.timestepping_order2)


		idstr += '_C'+str(self.timestep_size).zfill(6)

		if self.max_timesteps != -1:
			idstr += '_T'+str(self.max_timesteps).zfill(3)

		if self.rexi_method != '':
			if self.rexi_method == 'direct' or self.rexi_use_direct_solution:
				idstr += '_REXIDIR'
			else:
				if self.rexi_method == "file":
					idstr += '_REXIFIL'
					idstr += '_n'+str(self.rexi_file_n).zfill(8)
					idstr += '_h'+str(self.rexi_file_h)
					idstr += '_teabs'+str(self.rexi_file_test_abs).zfill(3)
					idstr += '_mxer'+str(self.rexi_file_max_error)

				elif self.rexi_method == "terry":
					idstr += '_REXITER'
					idstr += '_m'+str(self.rexi_m).zfill(8)
					idstr += '_h'+str(self.rexi_h)

				elif self.rexi_method == "ci":
					idstr += '_REXICI'
					idstr += '_n'+str(self.rexi_ci_n).zfill(8)
					if self.rexi_ci_max_real > 0:
						idstr += '_mr'+str(float(self.rexi_ci_max_real))
						idstr += '_mi'+str(float(self.rexi_ci_max_imag))
					else:
						idstr += '_sx'+str(float(self.rexi_ci_sx))
						idstr += '_sy'+str(float(self.rexi_ci_sy))
						idstr += '_mu'+str(float(self.rexi_ci_mu))
					idstr += '_pr'+str(self.rexi_ci_primitive)
					idstr += '_gfs'+str( "{:.4E}".format(self.rexi_ci_gaussian_filter_scale))
					idstr += '_gfd'+str( "{:.4E}".format(self.rexi_ci_gaussian_filter_dt_norm))
					idstr += '_gfe'+str( "{:.4E}".format(self.rexi_ci_gaussian_filter_exp_N))

				idstr += '_nrm'+str(self.rexi_normalization)
				idstr += '_hlf'+str(self.rexi_half_poles)
				idstr += '_bf'+str(self.rexi_beta_cutoff)
	
				#if self.plane_or_sphere == 'sphere':
				#idstr += '_pre'+str(self.rexi_sphere_preallocation)
				idstr += '_ext'+str(self.rexi_extended_modes).zfill(2)

				#idstr += '_rexithreadpar'+str(1 if self.rexi_thread_par else 0)

		if self.polvani_rossby >= 0:
			idstr += '_PR'+str(self.polvani_rossby)

		if self.polvani_froude >= 0:
			idstr += '_PF'+str(self.polvani_froude)

		if self.mode_res != -1:
			idstr += '_M'+str(self.mode_res).zfill(4)

		if self.phys_res != -1:
			idstr += '_N'+str(self.phys_res).zfill(4)

		return idstr
		

	def getRuntimeOptions(self):
		retval = ''
		retval += ' -g '+str(self.g)
		retval += ' -H '+str(self.h)
		retval += ' -f '+str(self.f)
		retval += ' -F '+str(self.f_sphere)
		retval += ' -a '+str(self.r)
		if self.mode_res != -1:
			retval += ' -M '+str(self.mode_res)
		if self.phys_res != -1:
			retval += ' -N '+str(self.phys_res)

		retval += ' --pde-id '+str(self.pde_id)
		retval += ' --staggering='+str(self.staggering)
		retval += ' -S '+str(self.spectralderiv)

		retval += ' -X '+str(self.domain_size)
		retval += ' -Y '+str(self.domain_size)
		
		retval += ' -s '+str(self.bench_id)
		retval += ' --benchmark='+str(self.benchmark_name)

		retval += ' -v '+str(self.verbosity)

		retval += ' --dt='+str(self.timestep_size)

		if self.max_timesteps != -1:
			retval += ' -T '+str(self.max_timesteps)

		retval += ' -o '+str(self.output_timestep_size)

		if self.output_timestep_size < 0:
			retval += ' -O -'	# deactivate file output

		if self.output_filename != '':
			retval += ' -O '+self.output_filename

		retval += ' -u '+str(self.viscosity)
		if self.viscosity_order > 0:
			retval += ' -U '+str(self.viscosity_order)
		retval += ' -t '+str(self.simtime)

		retval += ' --stability-checks='+str(self.stability_checks)

		if self.floating_point_output_digits >= 0:
			retval += ' -d '+str(self.floating_point_output_digits)

		retval += ' --use-linear-div='+str(self.uselineardiv)

		retval += ' --timestepping-method='+self.timestepping_method
		retval += ' --timestepping-order='+str(self.timestepping_order)
		retval += ' --timestepping-order2='+str(self.timestepping_order2)

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
					retval += ' --rexi-m='+str(self.rexi_m)
					retval += ' --rexi-h='+str(self.rexi_h)

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
					if self.rexi_ci_max_real > 0:
						retval += ' --rexi-ci-max-real='+str(self.rexi_ci_max_real)
						retval += ' --rexi-ci-max-imag='+str(self.rexi_ci_max_imag)
					else:
						retval += ' --rexi-ci-sx='+str(self.rexi_ci_sx)
						retval += ' --rexi-ci-sy='+str(self.rexi_ci_sy)
						retval += ' --rexi-ci-mu='+str(self.rexi_ci_mu)
					retval += ' --rexi-ci-primitive='+str(self.rexi_ci_primitive)
					retval += ' --rexi-ci-gaussian-filter-scale='+str(self.rexi_ci_gaussian_filter_scale)
					retval += ' --rexi-ci-gaussian-filter-dt-norm='+str(self.rexi_ci_gaussian_filter_dt_norm)
					retval += ' --rexi-ci-gaussian-filter-exp-N='+str(self.rexi_ci_gaussian_filter_exp_N)


		retval += ' --polvani-rossby='+str(self.polvani_rossby)
		retval += ' --polvani-froude='+str(self.polvani_froude)

		retval += ' --use-robert-functions='+str(self.use_robert_functions)

		retval += ' --compute-error='+str(self.compute_error)

		return retval

