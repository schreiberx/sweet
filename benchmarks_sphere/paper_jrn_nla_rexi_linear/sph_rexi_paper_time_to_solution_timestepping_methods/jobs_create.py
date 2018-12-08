#! /usr/bin/env python3

import os
import sys
import stat


class default_params:
	mode_res = 256

	output_timestep_size = -1

	f_sphere = 0

	timestepping_method = 'l_erk'
	timestepping_order = 1

	verbosity=0

	timestep_size = 0.001
	max_timesteps = 1

	normal_mode_analysis = 0

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

	simtime = 0.001 #math.inf

	rexi_par = 0
	postprocessing = 0

	compute_error = 0


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

make clean || exit 1
"""
#		content += """
#SCONS="scons --program=swe_sphere --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=release """+("--threading=off --rexi-thread-parallel-sum=enable" if p.rexi_par else "disable")+'"'+"""

		content += """
SCONS="scons --program=swe_sphere --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=release --threading=off"
echo "$SCONS"
$SCONS || exit 1
"""

		content += """
cd "$BASEDIR"
"""

		content += 'EXEC="$SWEETROOT/build/swe_sphere_spherespectral_spheredealiasing_libfft_gnu_release'
		content += ' -g '+str(self.g)
		content += ' -H '+str(self.h)
		content += ' -f '+str(self.f)
		content += ' -F '+str(self.f_sphere)
		content += ' -a '+str(self.r)
		content += ' -M '+str(self.space_res_spectral)

		content += ' --pde-id '+str(self.pde_id)

		content += ' -s '+str(self.bench_id)
		content += ' -v '+str(self.verbosity)

		content += ' -C '+str(-self.timestep_size)
		content += ' -T '+str(self.max_timesteps)

		content += ' -o '+str(self.output_timestep_size)
		content += ' -O -'	# deactivate file output
		content += ' -u '+str(self.viscosity)
		content += ' -t '+str(self.simtime)
		content += ' --nonlinear='+str(self.nonlinear)

		content += ' --timestepping-method='+self.timestepping_method
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
taskset -c 1 $EXEC || exit 1
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

		idstr += '_tsm_'+str(self.timestepping_method)
		idstr += '_tso'+str(self.timestepping_order)

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

p.use_robert_functions = 1
p.f_sphere = 0
p.rexi_extended_modes = 2



#
# Smaller values lead to no solution for the vort/div formulation
# See gaussian_ts_comparison_earth_scale_M16
# This shows that good results for RK can be computed with
# RK2 and a time step size of 200
#

p.timestep_size = 100
p.max_timesteps = 10
p.simtime = p.timestep_size*p.max_timesteps


for p.space_res_spectral in [64, 128, 256]:
	####################################
	# REXI
	####################################
	if True:
		p.timestepping_method = 'l_rexi'
		p.timestepping_order = 0

		p.rexi_m = 16

		p.gen_script('script'+p.create_job_id(), 'run.sh')

	p.rexi_m = 0

	####################################
	# RKn
	####################################
	for i in [1, 2, 3, 4]:
		p.timestepping_method = 'l_erk'
		p.timestepping_order = i

		p.gen_script('script'+p.create_job_id(), 'run.sh')

	####################################
	# LF2
	####################################
	if True:
		p.timestepping_method = 'l_lf'
		p.timestepping_order = 2

		p.gen_script('script'+p.create_job_id(), 'run.sh')

	####################################
	# IRK1
	####################################
	for i in [1]:
		p.timestepping_method = 'l_irk'
		p.timestepping_order = i

		p.gen_script('script'+p.create_job_id(), 'run.sh')

	####################################
	# CN
	####################################
	for i in [2]:
		p.timestepping_method = 'l_cn'
		p.timestepping_order = i

		p.gen_script('script'+p.create_job_id(), 'run.sh')


