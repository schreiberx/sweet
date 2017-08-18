#! /usr/bin/env python3

import os
import sys
import stat
import math


total_max_cores = 4096
max_cores_per_node = 16
total_max_nodes = total_max_cores/max_cores_per_node


class default_params:
	target_machine = ''
	#target_machine = 'yellowstone'

	prefix_string = ''

	plane_or_sphere = 'plane'

#	mode_res = 32
	mode_res = -1
	phys_res = -1

	output_timestep_size = 0.0001

	f_sphere = 1
	verbosity = 0

	timestepping_method = 'l_erk'
	timestepping_order = 1
	timestepping_order2 = 1

	timestep_size = 0.001
	max_timesteps = -1

	normal_mode_analysis = 0


	rexi_method = 'terry'
	rexi_file_n = 0
	rexi_file_h = 0
	rexi_file_test_abs = 0
	rexi_file_max_error = 0
	rexi_file_faf_dir = None

	rexi_m = 0
	rexi_l = 11
	rexi_h = 0.15

	rexi_half_poles = 1
	rexi_extended_modes = 0
	rexi_normalization = 0
	rexi_sphere_preallocation = 0
	rexi_use_direct_solution = 0


	def load_rexi_from_dict(self, d):
		if 'm' in d:
			self.rexi_m = d['m']

		if 'h' in d:
			self.rexi_h = d['h']

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

		if 'rexi_method' in d:
			self.rexi_method = d['rexi_method']

		if 'file_n' in d:
			self.rexi_file_n = d['file_n']

		if 'file_h' in d:
			self.rexi_file_h = d['file_h']

		if 'file_test_abs' in d:
			self.rexi_file_test_abs = d['file_test_abs']

		if 'file_max_error' in d:
			self.rexi_file_max_error = d['file_max_error']


	g = 1	# gravity
	h = 100000	# avg height

	f = 0.000072921	# \Omega coriolis effect
	#f = 2.0*f	# Coriolis effect for constant f (not multiplied with 2.0)

	r = 6371220	# radius
	domain_size = 1	# radius

	# 3: gaussian breaking dam
	# 4: geostrophic balance test case
	bench_id = 4
	use_robert_functions = 1

	pde_id = 0

	staggering = 0
	spectralderiv = 1
	nonlinear = 0
	viscosity = 0

	simtime = 0.001 #math.inf

	rexi_par = 0
	postprocessing = 0

	compute_error = 1

	# OpenMP threads in space
	par_space_threads = 1

	# MPI threads i time
	par_mpi_time_threads = 1


	def get_num_rexi_terms(self):
		rexi_terms = self.rexi_m+self.rexi_l

		if self.rexi_half_poles:
			rexi_terms = rexi_terms+1
		else:
			rexi_terms = 2*rexi_terms+1

		return rexi_terms


	def get_mpi_ranks_total(self):
		rexi_terms = self.get_num_rexi_terms()

#		if self.par_mpi_time_threads != 1:
#			if max_cores_per_node % self.par_space_threads != 0:
#				raise ValueError('Number of cores on node not evenly dividable by space threads')

		real_time_threads = min(rexi_terms, self.par_mpi_time_threads)

		# total number of used MPI ranks
		total_cores = self.par_space_threads*self.par_mpi_time_threads
		mpi_ranks_total = self.par_mpi_time_threads
		mpi_ranks_per_node = math.floor(max_cores_per_node/self.par_space_threads)

		return mpi_ranks_total


	def create_job_script(self, dirname):
		job_id = 'sweet_swe_sph_and_rexi'+self.create_job_id()

		rexi_terms = self.get_num_rexi_terms()

		if self.par_mpi_time_threads != 1:
			if max_cores_per_node % self.par_space_threads != 0:
				raise ValueError('Number of cores on node not evenly dividable by space threads')

		real_time_threads = min(rexi_terms, self.par_mpi_time_threads)

		# total number of used MPI ranks
		total_cores = self.par_space_threads*self.par_mpi_time_threads
		mpi_ranks_total = self.par_mpi_time_threads
		mpi_ranks_per_node = math.floor(max_cores_per_node/self.par_space_threads)

		program_bin_name=""
		program_bin_name=job_id

		sweetdir=os.path.normpath(os.getcwd()+'/../../')

		content = "#!/bin/bash\n"

		#
		# YELLOWSTONE:
		# Each node has 16 cores
		# 8 cores per socket
		# hyperthreading enabled
		#

		if self.target_machine == '':
			content += "\n"
		elif self.target_machine == 'yellowstone':
			content += """
#
# LSF batch script to run an MPI application
#
# YELLOW STONE SPECIFIC!!!
# https://www2.cisl.ucar.edu/resources/computational-systems/yellowstone/
#
#BSUB -P NCIS0002	# project code
#BSUB -W 02:00		# wall-clock time (hrs:mins)
#
#BSUB -n """+str(mpi_ranks_total)+"""	 number of tasks in job
#BSUB -R "span[ptile=16]"    # run 16 MPI tasks per node
#
#BSUB -outdir """+dirname+"""
#BSUB -J """+job_id+"""	# job name
#BSUB -o """+dirname+""".out  # output file name in which %J is replaced by the job ID
#BSUB -e """+dirname+""".out  # error file name in which %J is replaced by the job ID
#
## https://www2.cisl.ucar.edu/resources/computational-systems/yellowstone/using-computing-resources/queues-and-charges
#BSUB -q small
#

#
# More example job scripts:
# https://www2.cisl.ucar.edu/resources/computational-systems/yellowstone/using-computing-resources/running-jobs/platform-lsf-job-script-examples
#
""" # end yellowstone script
		else:
			print("Target machine "+str(self.target_machine)+" not supported")
			sys.exit(1)


		content += """

cd \""""+dirname+"""\"

BASEDIR="`pwd`"
rm -f ./prog_h_*
rm -f ./prog_u_*
rm -f ./prog_v_*

SWEETROOT=\""""+dirname+"""/../../../"
cd "$SWEETROOT"

pwd

# Always load local software
source ./local_software/env_vars.sh || exit 1

#make clean || exit 1

"""
		if self.target_machine == '':
			content += """
SCONS="scons --program=swe_plane_rexi --gui=disable --plane-spectral-space=enable --mode=release """+"--threading="+ ('omp' if not p.rexi_par else 'off') +" --rexi-thread-parallel-sum=" +('enable' if p.rexi_par else 'disable')+' -j 4"'+"""
echo "$SCONS"
$SCONS || exit 1

"""
		elif self.target_machine == 'yellowstone':
			pass

		else:
			print("Target machine "+str(self.target_machine)+" not supported")
			sys.exit(1)

		content += """
cd "$BASEDIR"
"""

		if self.rexi_par:
			content += 'EXEC="$SWEETROOT/build/swe_plane_rexi_planespectral_planedealiasing_rexipar_libfft_gnu_release'
		elif self.spectralderiv == 0:
			content += 'EXEC="$SWEETROOT/build/swe_plane_rexi_omp_libfft_gnu_release'
		else:
			content += 'EXEC="$SWEETROOT/build/swe_plane_rexi_planespectral_planedealiasing_omp_libfft_gnu_release'

		content += ' -g '+str(self.g)
		content += ' -H '+str(self.h)
		content += ' -f '+str(self.f)
		content += ' -F '+str(self.f_sphere)
		content += ' -a '+str(self.r)
		if self.mode_res != -1:
			content += ' -M '+str(self.mode_res)
		if self.phys_res != -1:
			content += ' -N '+str(self.phys_res)

		content += ' --pde-id '+str(self.pde_id)
		content += ' --staggering='+str(self.staggering)
		content += ' -S '+str(self.spectralderiv)

		content += ' -X '+str(self.domain_size)
		content += ' -s '+str(self.bench_id)

		content += ' -v '+str(self.verbosity)

		content += ' -C '+str(-self.timestep_size)

		if self.max_timesteps != -1:
			content += ' -T '+str(self.max_timesteps)

		content += ' -o '+str(self.output_timestep_size)

		if self.output_timestep_size < 0:
			content += ' -O -'	# deactivate file output

		content += ' -u '+str(self.viscosity)
		content += ' -t '+str(self.simtime)
		content += ' --nonlinear='+str(self.nonlinear)

		content += ' --timestepping-method='+self.timestepping_method
		content += ' --timestepping-order='+str(self.timestepping_order)
		content += ' --timestepping-order2='+str(self.timestepping_order2)

		content += ' --normal-mode-analysis-generation='+str(self.normal_mode_analysis)

		content += ' --rexi-method='+str(self.rexi_method)
		content += ' --rexi-half='+str(self.rexi_half_poles)
		content += ' --rexi-normalization='+str(self.rexi_normalization)
		content += ' --rexi-sphere-preallocation='+str(self.rexi_sphere_preallocation)
		content += ' --rexi-use-direct-solution='+str(self.rexi_use_direct_solution)
		content += ' --rexi-ext-modes='+str(self.rexi_extended_modes)

		# REXI Terry
		content += ' --rexi-m='+str(self.rexi_m)
		content += ' --rexi-h='+str(self.rexi_h)

		# REXI File
		content += ' --rexi-file-n='+str(self.rexi_file_n)
		content += ' --rexi-file-h='+str(self.rexi_file_h)
		content += ' --rexi-file-test-abs='+str(self.rexi_file_test_abs)
		content += ' --rexi-file-max-error='+str(self.rexi_file_max_error)
		if self.rexi_file_faf_dir != None:
			content += ' --rexi-file-faf-dir='+str(self.rexi_file_faf_dir)

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
		idstr = '_'+self.prefix_string


		idstr += '_bench'+str(self.bench_id)
#		idstr += '_nonlin'+str(self.nonlinear)

		idstr += '_g'+str(self.g)
		idstr += '_h'+str(self.h)
		idstr += '_f'+str(self.f)

		if self.plane_or_sphere == 'sphere':
			idstr += '_a'+str(self.r)

		idstr += '_u'+str(self.viscosity)

		idstr += '_pdeid'+str(self.pde_id)

		if self.plane_or_sphere == 'sphere':
			idstr += '_robert'+str(self.use_robert_functions)
			idstr += '_fsphere'+str(self.f_sphere)

#		idstr += '_t'+str(self.simtime).zfill(8)
#		idstr += '_o'+str(self.output_timestep_size).zfill(8)

		idstr += '_tsm_'+self.timestepping_method
		idstr += '_tso'+str(self.timestepping_order)
		idstr += '_tsob'+str(self.timestepping_order2)

		if True:
			if self.rexi_use_direct_solution:
				idstr += '_rexidirect'
			else:
				if self.rexi_method == "file":
					idstr += '_REXIFILE'
					idstr += '_filen'+str(self.rexi_file_n).zfill(8)
					idstr += '_fileh'+str(self.rexi_file_h)
					idstr += '_filetestabs'+str(self.rexi_file_test_abs).zfill(3)
					idstr += '_filemaxerr'+str(self.rexi_file_max_error)

				elif self.rexi_method == "terry":
					idstr += '_REXITERRY'
					idstr += '_rexim'+str(self.rexi_m).zfill(8)
					idstr += '_rexih'+str(self.rexi_h)

				elif self.rexi_method == "ci":
					idstr += '_REXICI'
					idstr += '_cim'+str(self.rexi_m).zfill(8)
					idstr += '_cih'+str(self.rexi_h)

				idstr += '_rexinorm'+str(self.rexi_normalization)
				idstr += '_rexihalf'+str(self.rexi_half_poles)
	
			if self.plane_or_sphere == 'sphere':
				idstr += '_rexiprealloc'+str(self.rexi_sphere_preallocation)
				idstr += '_rexiextmodes'+str(self.rexi_extended_modes).zfill(2)

#			idstr += '_rexipar'+str(1 if self.rexi_par else 0)

		idstr += '_C'+str(self.timestep_size).zfill(8)

		if self.max_timesteps != -1:
			idstr += '_Tn'+str(self.max_timesteps).zfill(3)

		if self.mode_res != -1:
			idstr += '_modes'+str(self.mode_res).zfill(4)

		if self.phys_res != -1:
			idstr += '_phys'+str(self.phys_res).zfill(4)

		return idstr


	def gen_script(self, dirname, scriptname):
		if not os.path.exists(dirname):
			os.makedirs(dirname)

		scriptname = 'run.sh'

		fullpath = dirname+'/'+scriptname
		print("WRITING "+fullpath)
		script_file = open(fullpath, 'w')
		script_file.write(self.create_job_script(os.getcwd()+'/'+dirname))
		script_file.close()

		st = os.stat(fullpath)
		os.chmod(fullpath, st.st_mode | stat.S_IEXEC)


p = default_params()


p.verbosity = 2

p.plane_or_sphere = 'plane'

p.mode_res = -1
p.phys_res = 128

p.bench_id = 1

p.rexi_sphere_preallocation = 0

p.g = 1
p.f = 1
p.h = 1
p.domain_size = 1

#p.viscosity = 0.0005
p.viscosity = 0.0

p.simtime = 0.1
p.output_timestep_size = p.simtime

timestep_size_reference = 0.0001
timestep_sizes = [0.0001*(2.0**i) for i in range(0, 11)]


# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']
#groups = ['ln2test']

#if len(sys.argv) < 5:
#	print("Usage: "+str(sys.argv[0])+" [group=l1/l2/ln1/ln2] [tsmethod] [order1] [order2] [use rexi direct solution]")
#	sys.exit(1)


if len(sys.argv) > 1:
	groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:
	# 1st order linear
	if group == 'l1':
		ts_methods = [
			['l_direct',	0,	0,	0],	# reference solution
			['l_erk',	1,	0,	0],
			['l_irk',	1,	0,	0],
			['l_rexi',	0,	0,	0],
		]

	# 2nd order linear
	if group == 'l2':
		ts_methods = [
			['l_direct',	0,	0,	0],	# reference solution
			['l_erk',	2,	0,	0],
			['l_cn',	2,	0,	0],
			['l_rexi',	0,	0,	0],
		]

	#	['lg_rexi_lc_erk_nt_sl_nd_erk',
	#	['l_rexi_ns_sl_nd_erk',

	# 1st order nonlinear
	if group == 'ln1':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			['l_erk_n_erk',		1,	1,	0],
			['l_irk_n_erk',		1,	1,	0],
			['ln_erk',		1,	1,	0],
			['l_rexi_n_erk',	1,	1,	0],
		]

	# 2nd order nonlinear
	if group == 'ln2':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			['l_cn_n_erk',		2,	2,	0],
			['l_erk_n_erk',		2,	2,	0],
			['l_irk_n_erk',		2,	2,	0],
			['ln_erk',		2,	2,	0],
			['l_rexi_n_erk',	2,	2,	0],
			['l_rexi_ns_sl_nd_erk',	2,	2,	0],
			['lg_rexi_lc_erk_nt_sl_nd_erk',	2,	2,	0],
		]


	# 4th order nonlinear
	if group == 'ln4':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			#['ln_etdrk',		4,	4,	1],	# reference solution

			['ln_etdrk',		4,	4,	1],
			['ln_erk',		4,	4,	0],
		]



	#
	# OVERRIDE TS methods
	#
	if len(sys.argv) > 4:
		ts_methods = [ts_methods[0]]+[[sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]]


	#
	# add prefix string to group benchmarks
	#
	prefix_string_template = group


	#
	# Reference solution
	#
	if True:
		tsm = ts_methods[0]

		p.prefix_string = prefix_string_template+'_ref'
		p.timestepping_method = tsm[0]
		p.timestepping_order = tsm[1]
		p.timestepping_order2 = tsm[2]
		p.rexi_use_direct_solution = tsm[3]

		if len(tsm) > 4:
			s = tsm[4]
			if 'timestep_size' in s:
				p.timestep_size = s['timestep_size']
		else:
			p.timestep_size = timestep_size_reference

		p.gen_script('script'+p.create_job_id(), 'run.sh')


	for tsm in ts_methods[1:]:
		for p.timestep_size in timestep_sizes:
			p.prefix_string = prefix_string_template

			p.timestepping_method = tsm[0]
			p.timestepping_order = tsm[1]
			p.timestepping_order2 = tsm[2]
			p.rexi_use_direct_solution = tsm[3]

			if len(tsm) > 4:
				s = tsm[4]
				p.load_rexi_from_dict(tsm[4])
				if 'timestep_size' in s:
					p.timestep_size = s['timestep_size']

			p.gen_script('script'+p.create_job_id(), 'run.sh')

