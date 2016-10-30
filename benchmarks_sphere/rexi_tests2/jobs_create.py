#! /usr/bin/env python

import os
import sys
import stat


class default_params:
	timestep_size = 0.01
	mode_res = 64
	output_timestep_size = 0.1
	rk_order = 4

	use_rexi = 1
	rexi_m = 256
	rexi_h = 0.2
	rexi_half_poles = 1
	rexi_extended_modes = 2

	rexi_use_coriolis_formulation = 1

	g = 1	# gravity
	h = 1	# avg height
	f = 0	# coriolis effect

	r = 10000	# radius


	bench_id = 3
	use_robert_functions = 1

	nonlinear = 0
	viscosity = 0

	simtime = 4

	rexi_par = 1
	postprocessing = 0


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
		content += """
SCONS="scons --program=swe_sph_and_rexi --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release """+("--threading=off --rexi-parallel-sum=enable" if p.rexi_par else "")+'"'+"""
echo "$SCONS"
$SCONS || exit 1
"""

		content += """
cd "$BASEDIR"
"""

		content += 'EXEC="$SWEETROOT/build/swe_sph_and_rexi_*_release'
		content += ' -g '+str(self.g)
		content += ' -H '+str(self.h)
		content += ' -f '+str(self.f)
		content += ' -s '+str(self.bench_id)
		content += ' -a '+str(self.r)
		content += ' -M '+str(self.mode_res)
		content += ' -C '+str(-self.timestep_size)
		content += ' -o '+str(self.output_timestep_size)
		content += ' -u '+str(self.viscosity)
		content += ' -t '+str(self.simtime)
		content += ' -R '+str(self.rk_order)
		content += ' --nonlinear='+str(self.nonlinear)
		content += ' --rexi='+str(self.use_rexi)
		content += ' --rexi-m='+str(self.rexi_m)
		content += ' --rexi-h='+str(self.rexi_h)
		content += ' --rexi-half='+str(self.rexi_half_poles)
		content += ' --use-robert-functions='+str(self.use_robert_functions)
		content += ' --rexi-ext-modes='+str(self.rexi_extended_modes)
		content += ' --rexi-use-coriolis-formulation='+str(self.rexi_use_coriolis_formulation)

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

		idstr += '_modes'+str(self.mode_res).zfill(3)
		idstr += '_bench'+str(self.bench_id)
		idstr += '_nonlin'+str(self.nonlinear)

		idstr += '_g'+str(self.g)
		idstr += '_h'+str(self.h)
		idstr += '_f'+str(self.f)
		idstr += '_a'+str(self.r)
		idstr += '_u'+str(self.viscosity)


		if self.use_rexi:
			idstr += '_rexim'+str(self.rexi_m).zfill(5)
			idstr += '_rexih'+str(self.rexi_h)
			idstr += '_rexihalf'+str(self.rexi_half_poles)
			idstr += '_rexiextmodes'+str(self.rexi_extended_modes).zfill(2)
			idstr += '_rexipar'+str(1 if self.rexi_par else 0)
			idstr += '_rexiusecoriolisformulation'+str(1 if self.rexi_use_coriolis_formulation else 0)
		else:
			idstr += '_RK'+str(self.rk_order)

		idstr += '_C'+str(self.timestep_size).zfill(8)
		idstr += '_t'+str(self.simtime).zfill(8)
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
p.simtime = 0.1
p.mode_res = 64

p.timestep_size = 10
p.output_timestep_size = 1000
p.simtime = 10000


####################################
####################################

if True:
	p.timestep_size = 10
	p.rexi_par = 1
	p.rexi_use_coriolis_formulation = 0

	for p.rexi_m in [64, 128, 256]:
		p.gen_script('script'+p.create_job_id(), 'run.sh')


if True:
	p.timestep_size = 100
	p.rexi_par = 1
	p.rexi_use_coriolis_formulation = 0

	for p.rexi_m in [64, 128, 256]:
		p.gen_script('script'+p.create_job_id(), 'run.sh')


####################################
# RK
####################################

if True:
	p.timestep_size = 10
	p.use_rexi = 0

	for p.rk_order in [1, 2, 3, 4]:
		p.gen_script('script'+p.create_job_id(), 'run.sh')


