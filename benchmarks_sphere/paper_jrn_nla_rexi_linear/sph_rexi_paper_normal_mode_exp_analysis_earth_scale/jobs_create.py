#! /usr/bin/env python3

import os
import sys
import stat
import math


class default_params:
#	mode_res = 32
	mode_res = 16

	output_timestep_size = 0.0001

	f_sphere = 1

	timestepping_method = 'l_erk'
	timestepping_order = 1

	timestep_size = 0.001
	max_timesteps = 1

	normal_mode_analysis = 1

	rexi_m = 256
	rexi_h = 0.15
	rexi_half_poles = 1
	rexi_extended_modes = 0
	rexi_normalization = 1
	rexi_sphere_preallocation = 1

	g = 1	# gravity
	h = 100000	# avg height

	f = 0.000072921	# \Omega coriolis effect
	#f = 2.0*f	# Coriolis effect for constant f (not multiplied with 2.0)

	r = 6371220	# radius

	# 3: gaussian breaking dam
	# 4: geostrophic balance test case
	bench_id = 4
	use_robert_functions = 1

	pde_id = 0

	nonlinear = 0
	viscosity = 0

	simtime = 0.001 #math.inf

	rexi_par = 1
	postprocessing = 0

	compute_error = 1


	def create_job_script(self):
		content = "#! /bin/bash\n"

		content += """
BASEDIR="`pwd`"
rm -f ./prog_h_*
rm -f ./prog_u_*
rm -f ./prog_v_*

SWEETROOT="../../../"
cd "$SWEETROOT"

# Always load local software
source ./local_software/env_vars.sh || exit 1

#make clean || exit 1
"""
		content += """
SCONS="scons --program=swe_sphere --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=release """+("--threading=off --rexi-thread-parallel-sum=enable" if p.rexi_par else "disable")+'"'+"""
echo "$SCONS"
$SCONS || exit 1
"""

		content += """
cd "$BASEDIR"
"""

		content += 'EXEC="$SWEETROOT/build/swe_sphere_spherespectral_spheredealiasing_rexipar_libfft_gnu_release'
		content += ' -g '+str(self.g)
		content += ' -H '+str(self.h)
		content += ' -f '+str(self.f)
		content += ' -F '+str(self.f_sphere)
		content += ' -a '+str(self.r)
		content += ' -M '+str(self.space_res_spectral)

		content += ' --pde-id '+str(self.pde_id)

		content += ' -s '+str(self.bench_id)

		content += ' -C '+str(-self.timestep_size)
		content += ' -T '+str(self.max_timesteps)

		content += ' -o '+str(self.output_timestep_size)
#		content += ' -O -'	# deactivate file output
		content += ' -u '+str(self.viscosity)
		content += ' -t '+str(self.simtime)
		content += ' --nonlinear='+str(self.nonlinear)

		content += ' --timestepping-method='+self.timestepping_method
		content += ' --timestepping-order='+str(self.timestepping_order)

		content += ' --normal-mode-analysis-generation='+str(self.normal_mode_analysis)

		content += ' --rexi-m='+str(self.rexi_m)
		content += ' --rexi-h='+str(self.rexi_h)
		content += ' --rexi-half='+str(self.rexi_half_poles)
		content += ' --rexi-normalization='+str(self.rexi_normalization)
		content += ' --rexi-sphere-preallocation='+str(self.rexi_sphere_preallocation)
		content += ' --rexi-ext-modes='+str(self.rexi_extended_modes)
		content += ' --use-robert-functions='+str(self.use_robert_functions)

		content += ' --compute-error='+str(self.compute_error)

		content += '"'

		content += "\n"
		content += """
echo "$EXEC"
$EXEC || exit 1
"""


		if self.postprocessing:
			content += """
../pp_plot_csv.py prog_h_*.csv
../pp_create_mp4.sh prog_h out_prog_h.mp4
"""

		return content


	def create_job_id(self):
		idstr = ''

#		idstr += '_modes'+str(self.space_res_spectral).zfill(3)
#		idstr += '_bench'+str(self.bench_id)
#		idstr += '_nonlin'+str(self.nonlinear)

		idstr += '_g'+str(self.g)
		idstr += '_h'+str(self.h)
		idstr += '_f'+str(self.f)
		idstr += '_a'+str(self.r)
		idstr += '_u'+str(self.viscosity)

		idstr += '_robert'+str(self.use_robert_functions)
		idstr += '_pdeid'+str(self.pde_id)
		idstr += '_fsphere'+str(self.f_sphere)

		idstr += '_C'+str(self.timestep_size).zfill(8)
		idstr += '_Tn'+str(self.max_timesteps).zfill(3)
		idstr += '_t'+str(self.simtime).zfill(8)
		idstr += '_o'+str(self.output_timestep_size).zfill(8)

		idstr += '_tsm_'+self.timestepping_method
		idstr += '_tso'+str(self.timestepping_order)

		if True:
			idstr += '_rexim'+str(self.rexi_m).zfill(8)
			idstr += '_rexih'+str(self.rexi_h)
			idstr += '_rexinorm'+str(self.rexi_normalization)
			idstr += '_rexihalf'+str(self.rexi_half_poles)
			idstr += '_rexiprealloc'+str(self.rexi_sphere_preallocation)
			idstr += '_rexiextmodes'+str(self.rexi_extended_modes).zfill(2)
#			idstr += '_rexipar'+str(1 if self.rexi_par else 0)

		return idstr


	def gen_script(self, dirname, scriptname):
		if not os.path.exists(dirname):
			os.makedirs(dirname)

		scriptname = 'run.sh'

		fullpath = dirname+'/'+scriptname
		print("WRITING "+fullpath)
		script_file = open(fullpath, 'w')
		script_file.write(self.create_job_script())
		script_file.close()

		st = os.stat(fullpath)
		os.chmod(fullpath, st.st_mode | stat.S_IEXEC)


p = default_params()



p.space_res_spectral = 16

# TODO: REPLACE THIS
#p.space_res_spectral = 12

p.normal_mode_analysis = 13

p.use_robert_functions = 1


default_timestep_size = 2000
default_timesteps = 1 #default_simtime/default_timestep_size

p.rexi_sphere_preallocation = 0


#for p.pde_id in [0, 1]:
for p.pde_id in [0]:

	for p.f_sphere in [0, 1]:

		if p.f_sphere == -1:
			p.sphere_rotating_coriolis_omega = 0
			p.f_sphere = 0

		elif p.f_sphere == 0:
			# f-sphere
			p.sphere_rotating_coriolis_omega = 0.000072921	# \Omega coriolis effect

		else:
			p.sphere_rotating_coriolis_omega = 0.000072921	# \Omega coriolis effect
			p.sphere_rotating_coriolis_omega = 2.0*p.f	# Constant f requires multiplication with 2.0


		####################################
		# REXI dt=defaut_timestep_size
		####################################
		for p.rexi_normalization in [1]:
			for p.rexi_extended_modes in [2]:
				p.timestepping_method = 'l_rexi'
				p.timestepping_order = 0

				p.timestep_size = default_timestep_size
				p.simtime = default_timestep_size*default_timesteps
				p.max_timesteps = default_timesteps

				for p.rexi_m in [2**i for i in range(4, 15)]:
					p.gen_script('script'+p.create_job_id(), 'run.sh')

		p.rexi_m = 0

		####################################
		# RK1
		####################################
		if 1:
			p.timestepping_method = 'l_erk'
			p.timestepping_order = 1

			p.timestep_size = default_timestep_size
			p.simtime = default_timestep_size*default_timesteps
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# RK2
		####################################
		if 1:
			p.timestepping_method = 'l_erk'
			p.timestepping_order = 2

			p.timestep_size = default_timestep_size
			p.simtime = default_timestep_size*default_timesteps
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# RK4
		####################################
		if 1:
			p.timestepping_method = 'l_erk'
			p.timestepping_order = 4

			p.timestep_size = default_timestep_size
			p.simtime = default_timestep_size*default_timesteps
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# IRK1
		####################################
		if 0:
			p.timestepping_method = 'l_irk'
			p.timestepping_order = 1

			p.timestep_size = default_timestep_size
			p.simtime = default_timestep_size*default_timesteps
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# IRK2
		####################################
		if 1:
			p.timestepping_method = 'l_cn'
			p.timestepping_order = 2

			p.timestep_size = default_timestep_size
			p.simtime = default_timestep_size*default_timesteps
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# LF2
		####################################
		if 0:
			p.timestepping_method = 'l_lf'
			p.timestepping_order = 2

			p.timestep_size = default_timestep_size
			p.simtime = default_timestep_size*default_timesteps
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


