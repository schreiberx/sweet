#! /usr/bin/env python

import os
import sys
import stat


class default_params:
#	mode_res = 32
	mode_res = 16

	output_timestep_size = 0.0001

	f_sphere = 1

	# 1: RK
	# 2: LF
	# 3: implicit RK, Crank-Nicolson
	timestepping_method = 1
	timestepping_order = 1

	timestep_size = 0.001
	max_timesteps = -1

	normal_mode_analysis = 1

	rexi_m = 256
	rexi_h = 0.15
	rexi_half_poles = 1
	rexi_extended_modes = 0

	g = 1	# gravity
	h = 100000	# avg height
	f = 0.00014584	# coriolis effect

	r = 6371220	# radius

	# 3: gaussian breaking dam
	# 4: geostrophic balance test case
	bench_id = 4
	use_robert_functions = 1

	pde_id = 0

	nonlinear = 0
	viscosity = 0

	simtime = -1

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
SCONS="scons --program=swe_sph_and_rexi --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=release """+("--threading=off --rexi-thread-parallel-sum=enable" if p.rexi_par else "disable")+'"'+"""
echo "$SCONS"
$SCONS || exit 1
"""

		content += """
cd "$BASEDIR"
"""

		#content += 'EXEC="$SWEETROOT/build/swe_sph_and_rexi_*_release'
		content += 'EXEC="$SWEETROOT/build/swe_sph_and_rexi_spherespectral_spheredealiasing_rexipar_libfft_gnu_release'
		content += ' -g '+str(self.g)
		content += ' -H '+str(self.h)
		content += ' -f '+str(self.f)
		content += ' -F '+str(self.f_sphere)
		content += ' -a '+str(self.r)
		content += ' -M '+str(self.mode_res)

		content += ' --pde-id '+str(self.pde_id)

		content += ' -s '+str(self.bench_id)

		content += ' -C '+str(-self.timestep_size)
		content += ' -T '+str(self.max_timesteps)

		content += ' -o '+str(self.output_timestep_size)
#		content += ' -O -'	# deactivate file output
		content += ' -u '+str(self.viscosity)
		content += ' -t '+str(self.simtime)
		content += ' --nonlinear='+str(self.nonlinear)

		content += ' --timestepping-method='+str(self.timestepping_method)
		content += ' --timestepping-order='+str(self.timestepping_order)

		content += ' --normal-mode-analysis-generation='+str(self.normal_mode_analysis)

		content += ' --rexi-m='+str(self.rexi_m)
		content += ' --rexi-h='+str(self.rexi_h)
		content += ' --rexi-half='+str(self.rexi_half_poles)
		content += ' --use-robert-functions='+str(self.use_robert_functions)
		content += ' --rexi-ext-modes='+str(self.rexi_extended_modes)

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

#		idstr += '_modes'+str(self.mode_res).zfill(3)
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

		idstr += '_tsm'+str(self.timestepping_method).zfill(3)
		idstr += '_tso'+str(self.timestepping_order)

#		if self.timestepping_method >= 100:
		if True:
			idstr += '_rexim'+str(self.rexi_m).zfill(8)
			idstr += '_rexih'+str(self.rexi_h)
			idstr += '_rexihalf'+str(self.rexi_half_poles)
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



p.mode_res = 16

# TODO: REPLACE THIS
#p.mode_res = 12

p.normal_mode_analysis = 13

p.use_robert_functions = 1


#
# Smaller values lead to no solution for the vort/div formulation
# See gaussian_ts_comparison_earth_scale_M16
# This shows that good results for RK can be computed with
# RK2 and a time step size of 200
#

default_timestep_size = 400
default_timesteps = 5 #default_simtime/default_timestep_size



#for p.pde_id in [0, 1]:
for p.pde_id in [1]:

	#for p.f_sphere in [-1, 0, 1]:
	for p.f_sphere in [0, 1]:

		if p.f_sphere == -1:
			p.f = 0
			p.f_sphere = 0
		else:
			p.f = 0.00014584


		####################################
		# REXI dt=defaut_timestep_size
		####################################
		for p.rexi_extended_modes in [2]:

			if 1:
				p.timestepping_method = 100
				p.timestepping_order = 0

				p.timestep_size = default_timestep_size
				p.max_timesteps = default_timesteps

				for p.rexi_m in [2**i for i in range(4, 8)]:
					p.gen_script('script'+p.create_job_id(), 'run.sh')

			if 1:
				p.timestepping_method = 100
				p.timestepping_order = 0

				p.timestep_size = default_timestep_size*default_timesteps
				p.max_timesteps = 1

				for p.rexi_m in [2**i for i in range(4, 11)]:
					p.gen_script('script'+p.create_job_id(), 'run.sh')

			p.rexi_m = 0


		####################################
		# RK4 tiny ts
		####################################
		if 0:
			p.timestepping_method = 1
			p.timestepping_order = 4

			p.timestep_size = default_timestep_size/16
			p.max_timesteps = default_timesteps*16

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# RK1
		####################################
		if 1:
			p.timestepping_method = 1
			p.timestepping_order = 1

			p.timestep_size = default_timestep_size
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# RK2
		####################################
		if 1:
			p.timestepping_method = 1
			p.timestepping_order = 2

			p.timestep_size = default_timestep_size
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# RK4
		####################################
		if 1:
			p.timestepping_method = 1
			p.timestepping_order = 4

			p.timestep_size = default_timestep_size
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# IRK1
		####################################
		if 0:
			p.timestepping_method = 3
			p.timestepping_order = 1

			p.timestep_size = default_timestep_size
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# IRK2
		####################################
		if 1:
			p.timestepping_method = 3
			p.timestepping_order = 2

			p.timestep_size = default_timestep_size
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


		####################################
		# LF2
		####################################
		if 0:
			p.timestepping_method = 2
			p.timestepping_order = 2

			p.timestep_size = default_timestep_size
			p.max_timesteps = default_timesteps

			p.gen_script('script'+p.create_job_id(), 'run.sh')


