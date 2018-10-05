#! /usr/bin/env python3

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
		self.user_script_exec_prefix = ''

		# Directly executed after the executable
		self.user_script_exec_suffix = ''

		#
		# REXI parallelization with threading?
		# True deactivates OpenMP space parallelization and uses OMP parallel for over the REXI terms
		#self.rexi_thread_par = False

		# Should the compile command be part of the job script
		self.compilecommand_in_jobscript = True

		# Filepath to jobscript
		self.p_job_scriptpath = None

		# Directory where to execute job
		self.p_job_dirpath = None

		# Filepath to write output of jobscript to
		self.p_job_stdout_filepath = None
		self.p_job_stderr_filepath = None

		self.p_job_picklepath = None

		# Accumulator for compile commands
		self.p_compilecommands_accum = []



	def print(self):
		self.print_info("compilecommand_in_jobscript: "+str(self.compilecommand_in_jobscript))
		self.print_info("p_job_scriptpath: "+str(self.p_job_scriptpath))
		self.print_info("p_job_dirpath: "+str(self.p_job_dirpath))
		self.print_info("p_job_stdout_filepath: "+str(self.p_job_stdout_filepath))
		self.print_info("p_job_stderr_filepath: "+str(self.p_job_stderr_filepath))
		self.print_info("p_job_picklepath: "+str(self.p_job_picklepath))



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

		if self.p_job_scriptpath == None:
			self.error("self.p_job_filepath == None")

		if self.p_job_dirpath == None:
			self.error("self.p_job_dirpath == None")

		if self.p_job_stdout_filepath == None:
			self.error("self.p_job_stdout_filepath == None")

		if self.p_job_stderr_filepath == None:
			self.error("self.p_job_stderr_filepath == None")



		self.compile.makeOptionsConsistent()

		content = ""

		# Setup job script generation
		self.platform_functions.jobscript_setup(self)

		# Now, assemble pieces together

		# Script header
		content += self.platform_functions.jobscript_get_header(self)
		content += self.user_script_header
		content += "# %SCRIPT_HEADER%\n"
			
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
				#print("INFO: Overriding environment variable "+i+"="+os.environ['SWEET_'+i])
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
cd \""""+self.p_job_dirpath+"""\"

"""

		# EXEC prefix
		content += self.platform_functions.jobscript_get_exec_prefix(self)
		content += self.user_script_exec_prefix
		content += "# %SCRIPT_EXEC_PREFIX%\n"

		# EXEC
		content += self.platform_functions.jobscript_get_exec_command(self)

		# EXEC suffix
		content += "# %SCRIPT_EXEC_SUFFIX%\n"
		content += self.user_script_exec_suffix
		content += self.platform_functions.jobscript_get_exec_suffix(self)

		content += self.user_script_footer
		content += self.platform_functions.jobscript_get_footer(self)
		content += "# %SCRIPT_FOOTER%\n"

		#
		# Add return and exit 0
		#
		# return in case that it's sourced into someone's bash
		# exit 0 is required in case that the last program (e.g. a copy) failed
		#
		content += "return 2>/dev/null; exit 0"

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

	def validate(self):
		if self.parallelization.max_wallclock_seconds != None and self.platform_resources.max_wallclock_seconds != None:
			if self.platform_resources.max_wallclock_seconds < self.parallelization.max_wallclock_seconds:
				self.info("platform_resources.max_wallclock_seconds: "+str(self.platform_resources.max_wallclock_seconds))
				self.info("parallelization.max_wallclock_seconds: "+str(self.parallelization.max_wallclock_seconds))
				self.error("Max. wallcock time exceeds platform's limit")


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



	def p_write_jobscript(self):

		# Validate before writing job script
		self.validate()

		# Generate script file content
		content = self.get_jobscript_content()

		# Write scriptfile
		script_file = open(self.p_job_scriptpath, 'w')
		self.info("Writing job file '"+self.p_job_scriptpath+"'")
		script_file.write(content)
		script_file.close()

		# Make script executable
		st = os.stat(self.p_job_scriptpath)
		os.chmod(self.p_job_scriptpath, st.st_mode | stat.S_IEXEC)



	def gen_jobscript_directory(self, i_job_dirpath):
		"""
		Generate a job script directory with
			run.sh:		Job script
			jobgeneration.pickle:	Information about job generation

		Parameters:
		--------------
		i_job_dirpath: str
			Directory to generate job script in.
			There is one directory per job!
		"""

		# Job directory
		self.p_job_dirpath = os.path.abspath(i_job_dirpath)

		# path to jobscript file
		self.p_job_scriptpath = self.p_job_dirpath+'/run.sh'

		# path to pickle file
		self.p_job_picklepath = self.p_job_dirpath+'/jobgeneration.pickle'

		# Determine output files
		self.p_job_stdout_filepath = self.p_job_dirpath+"/output.out"
		self.p_job_stderr_filepath = self.p_job_dirpath+"/output.err"


		# Generate directory of scriptfile if it does not exist
		if not os.path.exists(self.p_job_dirpath):
			os.makedirs(self.p_job_dirpath)

		self.p_write_jobscript()
		self.save_file(self.p_job_picklepath)



	def gen_script(self, dirname, scriptname):
		"""
		DEPRECATED
		DEPRECATED
		DEPRECATED
		"""
		self.info("WARNING: gen_script(...) is deprecated, use gen_jobscript_directory(...)")

		self.write_jobscript(dirname+'/'+scriptname)



	def write_jobscript(self, jobscript_filepath):
		"""
		DEPRECATED
		DEPRECATED
		DEPRECATED
		"""
		self.info("WARNING: write_jobscript(...) THIS FUNCTION IS DEPRECATED")

		# Job directory
		self.p_job_dirpath = os.path.abspath(os.path.dirname(jobscript_filepath))

		# path to jobscript file
		if os.path.basename(jobscript_filepath) != 'run.sh':
			raise Exception("Sorry, not supported anymore. Job script must be named 'run.sh'")

		self.p_job_scriptpath = self.p_job_dirpath+'/run.sh'

		# path to pickle file
		self.p_job_picklepath = self.p_job_dirpath+'/jobgeneration.pickle'

		# Determine output files
		self.p_job_stdout_filepath = self.p_job_dirpath+"/output.out"
		self.p_job_stderr_filepath = self.p_job_dirpath+"/output.err"

		# Generate directory of scriptfile if it does not exist
		if not os.path.exists(self.p_job_dirpath):
			os.makedirs(self.p_job_dirpath)

		self.p_write_jobscript()


	def getUniqueID(self, i_filter : list = []):
		"""
		Return a unique ID including *all* string and number attributes of this class

		i_filter:
			list of filter names to filter out from unique ID generation
		"""
		self.parallelization.dummy_setup_if_no_setup(self.platform_resources)
		unique_id = self.runtime.getUniqueID(self.compile, i_filter)

		s = self.compile.getUniqueParID(i_filter)
		if s != '':
			unique_id += '_'+s

		s = self.parallelization.getUniqueID(i_filter)
		if s != '':
			unique_id += '_'+s

		return unique_id



	def save_file(self, filename):
		import pickle

		def get_string_attributes(obj):
			attr_dict = {}
			obj_string_attributes = [a for a in dir(obj) if not a.startswith('__') and not callable(getattr(obj,a))]
			for attr in obj_string_attributes:
				a = getattr(obj, attr)
				if not isinstance(a, (float, int, str)):
					continue
				attr_dict[attr] = getattr(obj, attr)
			return attr_dict

		attr_dict = {}
		attr_dict['compile'] = get_string_attributes(self.compile)
		attr_dict['runtime'] = get_string_attributes(self.runtime)
		attr_dict['parallelization'] = get_string_attributes(self.parallelization)
		attr_dict['platforms_platform'] = get_string_attributes(self.platforms.platform)
		attr_dict['platform_resources'] = get_string_attributes(self.platform_resources)

		with open(filename, 'wb') as f:
			# Pickle the 'data' dictionary using the highest protocol available.
			pickle.dump(attr_dict, f, pickle.HIGHEST_PROTOCOL)

	def load_file(self, filename):
		import pickle

		def set_string_attributes(obj, attr_dict):
			for key, attr in attr_dict.items():
				setattr(obj, key, attr)

		with open(filename, 'rb') as f:
			# Pickle the 'data' dictionary using the highest protocol available.
			data = pickle.load(f)

			set_string_attributes(self.compile, data['compile'])
			set_string_attributes(self.runtime, data['runtime'])
			set_string_attributes(self.parallelization, data['parallelization'])
			set_string_attributes(self.platforms.platform, data['platforms_platform'])
			set_string_attributes(self.platform_resources, data['platform_resources'])



if __name__ == "__main__":
	p = SWEETJobGeneration()

	#s = p.get_jobscript_content()
	#p.info(s)

	s = p.getUniqueID()
	p.info(s)


	print("simtime: "+str(p.runtime.simtime))
	p.runtime.simtime = 666
	print("simtime: "+str(p.runtime.simtime))

	print("save_file()")
	p.save_file("/tmp/test.pickle")

	print("simtime: "+str(p.runtime.simtime))
	p.runtime.simtime = 123
	print("simtime: "+str(p.runtime.simtime))

	print("load_file()")
	p.load_file("/tmp/test.pickle")
	print("simtime: "+str(p.runtime.simtime))

	p.info("FIN")
