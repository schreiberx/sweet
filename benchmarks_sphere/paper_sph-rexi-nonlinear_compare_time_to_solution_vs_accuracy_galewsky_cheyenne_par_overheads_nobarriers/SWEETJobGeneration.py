#! /usr/bin/env python3

import os
import sys
import stat
import math
from SWEETCompileOptions import *
from SWEETRuntimeOptions import *
from SWEETClusterOptions import *


class SWEETJobGeneration:


	def __init__(self):
		self.compile = SWEETCompileOptions()
		self.runtime = SWEETRuntimeOptions()
		self.cluster = SWEETClusterOptions()

		self.user_script_header = ''
		self.user_script_footer = ''

		# Directly executed before the executable
		self.user_script_preprocess = ''

		# Directly executed after the executable
		self.user_script_postprocess = ''

		#
		# REXI parallelization with threading?
		# True deactivates OpenMP space parallelization and uses OMP parallel for over the REXI terms
		#self.rexi_thread_par = False




	def create_job_script(
			self,
			dirpath,	# directory path to generate benchmarks in
			dirname		# name of directory to store benchmark in
	):
		self.compile.makeOptionsConsistent()

		job_id = 'sweet_'+self.runtime.getUniqueID(self.compile)

		content, mpiexec_prefix = self.cluster.getScriptHeader('script'+self.runtime.getUniqueID(self.compile), self.runtime, self.compile, dirname)

		content += self.user_script_header
		content += """

cd \""""+dirpath+"""\"

BASEDIR="`pwd`"

SWEETROOT=\""""+dirpath+"""/../../../\"
cd "$SWEETROOT"
pwd

# Always load local software
#is this really the root?
if test -e ./local_software/env_vars.sh ; then
	source ./local_software/env_vars.sh || exit 1
else
	echo "Warning: changing SWEETROOT directory"	
	cd ..
	SWEETROOT="`pwd`"
	pwd
	source ./local_software/env_vars.sh || exit 1
fi

###make clean || exit 1
"""



		#
		# Setup compile options
		#

		if self.cluster.target_machine == '':
			content += """
SCONS="scons """+self.compile.getSConsParams()+' -j 4"'+"""
echo "$SCONS"
$SCONS || exit 1

"""

		elif self.cluster.target_machine == 'yellowstone':
			f = open('compile_yellowstone.sh', 'w')
			f.write("#! /bin/bash\n")
			f.write("\n")
			f.write("scons "+self.compile.getSConsParams()+'\n')
			f.write("\n")

		elif self.cluster.target_machine in ['cheyenne', 'cheyenne_impi', 'cheyenne_openmpi', 'mac-login-amd', 'mac-login-intel']:
			fn = 'compile_'+self.cluster.target_machine+'.sh'
			f = open(fn, 'w')
			f.write("#! /bin/bash\n")
			f.write("\n")
			f.write("SWEETROOT=\""+dirpath+"/../../../\"\n")
			f.write("cd \"$SWEETROOT\"\n")
			f.write("\n")

			if self.cluster.target_machine == 'cheyenne_impi':
				# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/intel-mpi-and-open-mpi
				f.write("module load impi\n")
				f.write("\n")

			if self.cluster.target_machine == 'cheyenne_openmpi':
				# https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/intel-mpi-and-open-mpi
				f.write("module load openmpi\n")
				f.write("\n")

			if self.compile.mkl == 'enable':
				f.write("module load mkl\n")
				f.write("\n")

			f.write("scons "+self.compile.getSConsParams()+'\n')
			f.write("\n")
			os.chmod(fn, 0o755)

			#print("COMPILE WITH: scons "+self.compile.getSConsParams()+' -j 4')

		else:
			content += "\n"
			content += "scons "+self.compile.getSConsParams()+"\n"
			content += "\n"

		content += """

cd "$BASEDIR"
pwd

"""

		content += self.user_script_preprocess

		content += 'EXEC="$SWEETROOT/build/'+self.compile.getProgramName()+' '
		content += self.runtime.getRuntimeOptions()
		content += '"'

		content += "\n"
		content += """

echo "$EXEC"
pwd
#ln -s "$SWEETROOT/data/" "$BASEDIR/data"   #Symlink for GUI directory, if necessary
"""+mpiexec_prefix+"""$EXEC || exit 1

"""
		content += self.user_script_postprocess

		content += self.user_script_footer

		return content



	def gen_script(self, dirname, scriptname):
		if not os.path.exists(dirname):
			os.makedirs(dirname)

		scriptname = 'run.sh'

		fullpath = dirname+'/'+scriptname
		print("WRITING "+fullpath)
		script_file = open(fullpath, 'w')
		script_file.write(self.create_job_script(os.getcwd()+'/'+dirname, dirname))
		script_file.close()

		st = os.stat(fullpath)
		os.chmod(fullpath, st.st_mode | stat.S_IEXEC)

