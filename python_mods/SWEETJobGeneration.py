
import os
import sys
import stat

from SWEETCompileOptions import *
from SWEETRuntimeOptions import *
from SWEETPlatforms import *
from InfoError import *
from SWEETParallelization import *

__all__ = ['SWEETJobGeneration']

class SWEETJobGeneration(InfoError):


	def __init__(self, platform_id_override = None):
		InfoError.__init__(self, "SWEETJobGeneration")

		self.sweetroot = os.environ.get('SWEET_ROOT')
		if self.sweetroot == None:
			self.error("SWEET environment variables not loaded!")


		# Setup all options
		self.compile : SWEETCompileOptions = SWEETCompileOptions()
		self.runtime : SWEETRuntimeOptions = SWEETRuntimeOptions()
		self.parallelization : SWEETParallelization = SWEETParallelization()

		self.platforms = SWEETPlatforms(platform_id_override)
		self.platform_functions = self.platforms.functions

		self.platform_resources = self.platform_functions.get_platform_resources()

		# Setup all hardware print_information which also ensures consistency
		self.platform_resources.setup()
		self.platform_resources.print()

		# Userdefined data for script header
		self.user_script_header = ''

		# Userdefined data for script footer
		self.user_script_footer = ''

		# Directly executed before the executable
		self.user_script_preprocess = ''

		# Directly executed after the executable
		self.user_script_postprocess = ''

		#
		# REXI parallelization with threading?
		# True deactivates OpenMP space parallelization and uses OMP parallel for over the REXI terms
		#self.rexi_thread_par = False

		# Should the compile command be part of the job script
		self.compilecommand_in_jobscript = True

		# Filepath to jobscript
		self.p_jobscript_filepath = None

		# Directory where to execute job
		self.p_jobscript_dirpath = None

		# Filepath to write output of jobscript to
		self.p_jobscript_stdout_filepath = None
		self.p_jobscript_stderr_filepath = None

		# Accumulator for compile commands
		self.p_compilecommands_accum = []



	def print(self):
		self.print_info("compilecommand_in_jobscript: "+str(self.compilecommand_in_jobscript))
		self.print_info("p_jobscript_filepath: "+str(self.p_jobscript_filepath))
		self.print_info("p_jobscript_dirpath: "+str(self.p_jobscript_dirpath))
		self.print_info("p_jobscript_stdout_filepath: "+str(self.p_jobscript_stdout_filepath))
		self.print_info("p_jobscript_stderr_filepath: "+str(self.p_jobscript_stderr_filepath))



	def setup_parallelization(self, parallelization_dim_list):
		"""
		Setup the parallelization configuration

		Parameters
		----------
		parallelization_dim_list: list
			list with SWEETParallelizationDimOptions with parallelization
			print_information along each dimension
		"""
		self.parallelization.setup(parallelization_dim_list, self.platform_resources)


	def get_program_exec(self):
		"""
		Return program executable and options

		This is accessed by the platform specific scripts to get the executable
		"""
		return self.sweetroot+"/build/"""+self.compile.getProgramName()+""" """+self.runtime.getRuntimeOptions()



	def get_jobscript_content(self):
		"""
		Generate the job script based on configuration
		"""

		# Run Dummy setup in case that no setup was done so far!
		self.parallelization.dummy_setup_if_no_setup(self.platform_resources)

		if self.p_jobscript_filepath == None:
			self.error("self.p_jobscript_filepath == None")

		if self.p_jobscript_dirpath == None:
			self.error("self.p_jobscript_dirpath == None")

		if self.p_jobscript_stdout_filepath == None:
			self.error("self.p_jobscript_stdout_filepath == None")

		if self.p_jobscript_stderr_filepath == None:
			self.error("self.p_jobscript_stderr_filepath == None")


		self.compile.makeOptionsConsistent()

		content = ""

		# Setup job script generation
		self.platform_functions.jobscript_setup(self)

		# Now, assemble pieces together

		# Script header
		content += self.platform_functions.jobscript_get_header(self)
		content += self.user_script_header
			
		content += """

# Provide platform ID helper for scripts.
# This makes platform detection on the compute nodes easier since they might have mainly numerical identifiers
export SWEET_PLATFORM=\""""+self.platforms.platform_id+"""\"

# Loading SWEET environment variables
cd \""""+self.sweetroot+"""\"
source ./local_software/env_vars.sh \""""+os.path.normpath(self.platforms.platform_dirpath+"/env_vars.sh")+"""\" || exit 1

"""
		#
		# Override compilers in environment variables
		#
		override_list = ['CC', 'CXX', 'F90', 'MPICC', 'MPICXX', 'MPIF90']

		for i in override_list:
			if 'SWEET_'+i in os.environ:
				print("INFO: Overriding environment variable "+i+"="+os.environ['SWEET_'+i])
				content += "export "+i+"="+os.environ['SWEET_'+i]+"\n"

		# Compile in main SWEET directory if requested
		c = self.platform_functions.jobscript_get_compile_command(self)
		if self.compilecommand_in_jobscript:
			content += c
		else:
			if c not in self.p_compilecommands_accum:
				self.p_compilecommands_accum.append(c)


		# Change to execution directory
		content += """

# chdir to execution directory
cd \""""+self.p_jobscript_dirpath+"""\"

"""


		# EXEC prefix
		content += self.platform_functions.jobscript_get_exec_prefix(self)
		content += self.user_script_preprocess

		# EXEC
		content += self.platform_functions.jobscript_get_exec_command(self)

		# EXEC suffix
		content += self.user_script_postprocess
		content += self.platform_functions.jobscript_get_exec_suffix(self)

		content += self.user_script_footer
		content += self.platform_functions.jobscript_get_footer(self)

		return content


	def get_compilecommands_content(self):
		"""
		Return the content of a scriptfile to compile all programs for which the JobGeneration script was used for
		"""
		content = """#! /bin/bash

# Loading SWEET environment variables
cd \""""+self.sweetroot+"""\"
source ./local_software/env_vars.sh \""""+os.path.normpath(self.platforms.platform_dirpath+"/env_vars.sh")+"""\" || exit 1

"""
		content += "\n".join(self.p_compilecommands_accum)

		# Empty accumulator
		self.p_compilecommands_accum = []
		return content



	def write_compilecommands(self, compilecommands_filename):
		compilecommands_filepath = os.path.abspath(compilecommands_filename)

		# Create directory for job script if it doesn't yet exist
		compilecommands_dirpath = os.path.dirname(compilecommands_filepath)

		# Generate script file content
		content = self.get_compilecommands_content()

		self.info("Compile commands file '"+compilecommands_filepath+"'")

		# Generate directory of scriptfile is not exist
		if not os.path.exists(compilecommands_dirpath):
			os.makedirs(compilecommands_dirpath)

		# Write scriptfile
		script_file = open(compilecommands_filepath, 'w')
		script_file.write(content)
		script_file.close()

		# Set permissions
		st = os.stat(compilecommands_filepath)
		os.chmod(compilecommands_filepath, st.st_mode | stat.S_IEXEC)



	def write_jobscript(self, jobscript_filepath):

		if self.parallelization.max_wallclock_seconds != None and self.platform_resources.max_wallclock_seconds != None:
			if self.platform_resources.max_wallclock_seconds < self.parallelization.max_wallclock_seconds:
				self.info("platform_resources.max_wallclock_seconds: "+str(self.platform_resources.max_wallclock_seconds))
				self.info("parallelization.max_wallclock_seconds: "+str(self.parallelization.max_wallclock_seconds))
				self.error("Max. wallcock time exceeds platform's limit")

		if jobscript_filepath[0] == '/':
			self.p_jobscript_filepath = jobscript_filepath
		else:
			self.p_jobscript_filepath = os.path.abspath(jobscript_filepath)

		# Create directory for job script if it doesn't yet exist
		self.p_jobscript_dirpath = os.path.dirname(self.p_jobscript_filepath)

		# Execute in directory where job is generated in

		# Setup default output
		self.p_jobscript_stdout_filepath = self.p_jobscript_dirpath+"/output.out"
		self.p_jobscript_stderr_filepath = self.p_jobscript_dirpath+"/output.err"

		# Generate script file content
		content = self.get_jobscript_content()

		self.info("Job file '"+self.p_jobscript_filepath+"'")

		# Generate directory of scriptfile is not exist
		if not os.path.exists(self.p_jobscript_dirpath):
			os.makedirs(self.p_jobscript_dirpath)

		# Write scriptfile
		script_file = open(self.p_jobscript_filepath, 'w')
		script_file.write(content)
		script_file.close()

		# Set permissions
		st = os.stat(self.p_jobscript_filepath)
		os.chmod(self.p_jobscript_filepath, st.st_mode | stat.S_IEXEC)




	def getUniqueID(self):
		self.parallelization.dummy_setup_if_no_setup(self.platform_resources)
		return self.runtime.getUniqueID(self.compile)+'_'+self.parallelization.getUniqueID()



	# Keep this for backward compatibility
	def gen_script(self, dirname, scriptname):
		#self.info("WARNING: gen_script(...) is deprecated, use write_jobscript(...)")

		self.write_jobscript(dirname+'/'+scriptname)



if __name__ == "__main__":
	p = SWEETJobGeneration()

	#s = p.get_jobscript_content()
	#p.info(s)

	s = p.getUniqueID()
	p.info(s)

	p.info("FIN")
