#! /usr/bin/env python3

import os
import sys
import stat


class default_params:
	timestep_size = 0.001
	mode_res = 64
	output_timestep_size = 0.01
	rk_order = 4

	timestepping_method = 1
	rexi_m = 256
	rexi_h = 0.15
	rexi_half_poles = 1
	sphere_extended_modes = 2
	rexi_normalization = 1

	pde_id = 0

	g = 1	# gravity
	h = 1	# avg height
	f = 1	# coriolis effect

	r = 1	# radius

	# 10: geostrophic balance test case
	bench_id = 10
	use_robert_functions = 1

	nonlinear = 0
	viscosity = 0

	simtime = 1

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
SCONS="scons --program=swe_sphere --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release """+("--threading=off --rexi-thread-parallel-sum=enable" if p.rexi_par else "")+' -j 4"'+"""
echo "$SCONS"
$SCONS || exit 1
"""

		content += """
cd "$BASEDIR"
"""

		#content += 'EXEC="$SWEETROOT/build/swe_sph_and_rexi_*_release'
		content += 'EXEC="$SWEETROOT/build/swe_sphere_spherespectral_spheredealiasing_rexipar_libfft_gcc_release'

		if self.gravitation>= 0:
			content += ' -g '+str(self.g)

		if self.h0 >= 0:
			content += ' -H '+str(self.h)

		if self.sphere_rotating_coriolis_omega >= 0:
			content += ' -f '+str(self.f)

		if self.sphere_radius >= 0:
			content += ' -a '+str(self.r)

		content += ' -s '+str(self.bench_id)
		content += ' -M '+str(self.space_res_spectral)
		content += ' -C '+str(-self.timestep_size)
		content += ' -o '+str(self.output_timestep_size)
		content += ' -O -'	# deactivate file output
		content += ' -u '+str(self.viscosity)
		content += ' -t '+str(self.max_simulation_time)
		content += ' -R '+str(self.rk_order)
		content += ' --pde-id '+str(self.pde_id)
		content += ' --nonlinear='+str(self.nonlinear)
		content += ' --timestepping-method='+str(self.timestepping_method)
		content += ' --timestepping-order='+str(self.timestepping_order)
		content += ' --rexi-m='+str(self.rexi_m)
		content += ' --rexi-h='+str(self.rexi_h)
		content += ' --rexi-half='+str(self.rexi_half_poles)
		content += ' --rexi-normalization='+str(self.rexi_normalization)
		content += ' --rexi-ext-modes='+str(self.sphere_extended_modes)
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

		idstr += '_modes'+str(self.space_res_spectral).zfill(3)
		idstr += '_bench'+str(self.bench_id)
		idstr += '_nonlin'+str(self.nonlinear)

		idstr += '_g'+str(self.g)
		idstr += '_h'+str(self.h)
		idstr += '_f'+str(self.f)
		idstr += '_a'+str(self.r)
		idstr += '_u'+str(self.viscosity)

		idstr += '_pdeid'+str(self.pde_id)

		idstr += '_tsm'+str(self.timestepping_method)
		idstr += '_tso'+str(self.timestepping_order)

#		if self.use_rexi:
		if True:
			idstr += '_rexim'+str(self.rexi_m).zfill(5)
			idstr += '_rexih'+str(self.rexi_h)
			idstr += '_rexinorm'+str(self.rexi_normalization)
			idstr += '_rexihalf'+str(self.rexi_half_poles)
			idstr += '_rexiextmodes'+str(self.sphere_extended_modes).zfill(2)
			idstr += '_rexipar'+str(1 if self.rexi_par else 0)

		idstr += '_C'+str(self.timestep_size).zfill(8)
		idstr += '_t'+str(self.max_simulation_time).zfill(8)
		idstr += '_o'+str(self.output_timestep_size).zfill(8)

		idstr += '_robert'+str(self.use_robert_functions)

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


####################################
# REXI
####################################

p.max_simulation_time = 0.1
p.bench_id = 10	# Geostrophic balance benchmark

p.output_timestep_size = 0.01

if False:
	self.gravitation= -1
	self.h0 = -1
	self.sphere_rotating_coriolis_omega = -1
	self.sphere_radius = -1

if True:
	# 10 times larger than RK4 time step size
	p.timestepping_method = 'l_rexi'
	p.timestepping_order = 1

	p.timestep_size = 0.01
	p.rexi_par = 1

	for p.rexi_half_poles in [0, 1]:
		for p.rexi_normalization in [1,0]:
			for p.sphere_extended_modes in [2]:
				for p.space_res_spectral in [64]:
					for p.rexi_m in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
						p.gen_script('script'+p.create_job_id(), 'run.sh')



####################################
# RK
####################################

if True:
	p.timestepping_method = 'l_erk'
	p.timestepping_order = 2

	p.timestep_size = 0.01

	for p.space_res_spectral in [64]:
		p.gen_script('script'+p.create_job_id(), 'run.sh')


