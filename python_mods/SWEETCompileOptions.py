#! /usr/bin/env python3

#
# This file describes all possible compile options
#
# The underlying concept is to separate compile options
# and the unique Id generation into this to make it reusable
# by other scripts.
#
# This also leads to a clean separation between compile and runtime options.
#

import sys
import os
from importlib import import_module
from InfoError import *

import subprocess


__all__ = ['SWEETCompileOptions']


def _exec_command(command):
	process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	# combine stdout and stderr
	out = out+err
	out = out.decode("utf-8")
	out = out.replace("\r", "")
	if len(out) > 0:
		if out[-1] in ["\n", "\r"]:
			out = out[:-1]
	return out


class SWEETCompileOptions(InfoError):


	def __init__(self):
		InfoError.__init__(self, "SWEETCompileOptions")

		# Program or unit test
		self.program = ""
		self.unit_test = ""

		# Compile options
		self.mode = 'release'
		self.compiler = 'gnu'

		self.debug_symbols = 'enable'
		self.simd = 'enable'
		self.mic = 'disable'

		self.gxx_toolchain = ''
		self.cxx_flags = ''
		self.ld_flags = ''

		self.fortran_source = 'disable'

		self.lapack = 'enable'

		self.program_binary_name = ''

		# LLVM overrides
		self.llvm_gnu_override = False
		self.llvm_omp_override = False


		# Parallelization
		self.sweet_mpi = 'disable'
		self.threading = 'omp'
		self.rexi_thread_parallel_sum = 'disable'

		self.rexi_timings = 'disable'

		# Additional barriers to overcome issues of turbo boost
		self.rexi_timings_additional_barriers = 'disable'

		# Memory allocator
		if _exec_command('uname -s') == "Darwin":
			# Deactivate efficient NUMA block allocation on MacOSX systems (missing numa.h file)
			self.numa_block_allocator = 0
		else:
			self.numa_block_allocator = 2

		# Program / Unit test
		self.program = ''
		self.unit_test = ''

		# PinT
		self.parareal = 'none'
		self.libpfasst = 'disable'
		self.pfasst_cpp = 'disable'

		# Eigen library
		self.eigen = 'disable'

		# Libraries
		self.libfft = 'disable'
		self.libsph = 'disable'
		self.mkl = 'disable'

		# Features
		self.plane_spectral_space = 'enable'
		self.plane_spectral_dealiasing = 'enable'
		self.sphere_spectral_space = 'disable'
		self.sphere_spectral_dealiasing = 'enable'
		self.libxml = 'disable'

		# GUI
		self.gui = 'disable'

		#self.quadmath = 'enable'
		self.quadmath = 'disable'
		pass


	def getSConsParams(self):
		retval = ''
		retval += ' --mode='+self.mode
		retval += ' --compiler='+self.compiler
		retval += ' --debug-symbols='+("enable" if self.debug_symbols else "disable")
		retval += ' --simd='+self.simd
		retval += ' --mic='+self.mic

		if self.gxx_toolchain != '':
			retval += ' --gxx-toolchain='+self.gxx_toolchain
		if self.cxx_flags != '':
			retval += ' --cxx-flags='+self.cxx_flags
		if self.ld_flags:
			retval += ' --ld-flags='+self.ld_flags

		retval += ' --fortran-source='+self.fortran_source

		retval += ' --lapack='+self.lapack

		retval += ' --program-binary-name='+self.program_binary_name

		# Parallelization
		retval += ' --sweet-mpi='+self.sweet_mpi
		retval += ' --threading='+self.threading
		retval += ' --rexi-thread-parallel-sum='+self.rexi_thread_parallel_sum
		retval += ' --rexi-timings='+self.rexi_timings
		retval += ' --rexi-timings-additional-barriers='+self.rexi_timings_additional_barriers

		# Memory allocator
		retval += ' --numa-block-allocator='+str(self.numa_block_allocator)

		# Program / Unit test
		if self.program != '':
			retval += ' --program='+self.program
		if self.unit_test != '':
			retval += ' --unit-test='+self.unit_test

		# PinT
		retval += ' --parareal='+self.parareal
		retval += ' --libpfasst='+self.libpfasst
		retval += ' --pfasst-cpp='+self.pfasst_cpp

		retval += ' --eigen='+self.eigen

		# Libraries
		retval += ' --libfft='+self.libfft
		retval += ' --libsph='+self.libsph
		retval += ' --mkl='+self.mkl

		# Features
		retval += ' --plane-spectral-space='+self.plane_spectral_space
		retval += ' --plane-spectral-dealiasing='+self.plane_spectral_dealiasing
		retval += ' --sphere-spectral-space='+self.sphere_spectral_space
		retval += ' --sphere-spectral-dealiasing='+self.sphere_spectral_dealiasing
		retval += ' --libxml='+self.libxml

		# GUI
		retval += ' --gui='+self.gui

		# Activate quadmath
		retval += ' --quadmath='+self.quadmath

		return retval



	def makeOptionsConsistent(self):

		if self.libpfasst == 'enable':
			self.fortran_source = 'enable'

		if self.plane_spectral_space == 'enable':
			self.libfft = 'enable'
		else:
			self.plane_spectral_dealiasing = 'disable'

		if self.sphere_spectral_space == 'enable':
			self.libsph = 'enable'
		else:
			self.sphere_spectral_dealiasing = 'disable'

		if self.libsph == 'enable':
			# activate linking with libfft!
			self.libfft = 'enable'

		return


	def sconsProcessOptions(self):
		scons = import_module("SCons.Script")

#		scons.AddOption(	'--xml-config',
#				dest='xml_config',
#				type='string',
#				default='',
#				help='xml configuration file with compile options'
#		)
#		env['xml_config'] = scons.GetOption('xml_config')
#
#		if env['xml_config'] != '':
#			self.env = CompileXMLOptions.load(env['xml_config'], env)



		scons.AddOption(      '--mode',
				dest='mode',
				type='choice',
				choices=['debug', 'release'],
				default='release',
				help='specify compiler to use: debug, release [default: %default]'
		)
		self.mode = scons.GetOption('mode')


		scons.AddOption(	'--compiler',
				dest='compiler',
				type='choice',
				choices=['gnu', 'intel', 'llvm', 'pgi'],
				default='gnu',
				help='specify compiler to use: gnu, intel, llvm, pgi [default: %default]'
		)
		self.compiler = scons.GetOption('compiler')

		scons.AddOption(	'--numa-block-allocator',
				dest='numa_block_allocator',
				type='choice',
				choices=['0', '1', '2', '3'],
				default=str(self.numa_block_allocator),
				help='Specify allocation method to use: 0: default system\'s malloc, 1: allocation with NUMA granularity, 2: allocation with thread granularity, 3: allocation with non-NUMA granularity [default: %default]'
		)
		self.numa_block_allocator = int(scons.GetOption('numa_block_allocator'))


		scons.AddOption(	'--gxx-toolchain',
				dest='gxx-toolchain',
				type='string',
				default='',
				help='specify gcc toolchain for intel and llvm compiler, e.g. g++-4.6, default: deactivated'
		)
		self.gxx_toolchain = scons.GetOption('gxx-toolchain')


		scons.AddOption(	'--simd',
				dest='simd',
				type='choice',
				choices=['enable', 'disable'],
				default='enable',
				help="Use SIMD for operations such as folding [default: %default]"
		)
		self.simd = scons.GetOption('simd')


		scons.AddOption(	'--eigen',
				dest='eigen',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Activate utilization of Eigen library [default: %default]"
		)
		self.eigen = scons.GetOption('eigen')



		scons.AddOption(	'--pfasst-cpp',
				dest='pfasst_cpp',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Activate utilization of PFASST++ (C++ version) [default: %default]"
		)
		self.pfasst_cpp = scons.GetOption('pfasst_cpp')


		scons.AddOption(	'--libpfasst',
				dest='libpfasst',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Activate utilization of libPFASST (FOortran version) [default: %default]"
		)
		self.libpfasst = scons.GetOption('libpfasst')



		scons.AddOption(	'--debug-symbols',
				dest='debug_symbols',
				type='choice',
				choices=['enable', 'disable'],
				default='enable',
				help="Create binary with debug symbols [default: %default]"
		)
		self.debug_symbols = scons.GetOption('debug_symbols')



		scons.AddOption(	'--plane-spectral-space',
				dest='plane_spectral_space',
				type='choice',
				choices=['enable', 'disable'],
				default='enable',
				help="Activate spectral space for data on the plane (2D FFT) [default: %default]"
		)
		self.plane_spectral_space = scons.GetOption('plane_spectral_space')

		scons.AddOption(	'--sphere-spectral-space',
				dest='sphere_spectral_space',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Activate spectral space for data on the sphere (Spherical Harmonics) [default: %default]"
		)
		self.sphere_spectral_space = scons.GetOption('sphere_spectral_space')


		scons.AddOption(	'--fortran-source',
				dest='fortran_source',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Activate linking with Fortran source [default: %default]"
		)
		self.fortran_source = scons.GetOption('fortran_source')


		scons.AddOption(	'--lapack',
				dest='lapack',
				type='choice',
				choices=['enable', 'disable'],
				default='enable',
				help="Enable lapack [default: %default]"
		)
		self.lapack = scons.GetOption('lapack')



		scons.AddOption(	'--libfft',
				dest='libfft',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Enable compiling and linking with FFT library [default: %default]"
		)
		self.libfft = scons.GetOption('libfft')


		scons.AddOption(	'--libsph',
				dest='libsph',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Enable compiling and linking with SPH library [default: %default]"
		)
		self.libsph = scons.GetOption('libsph')


		scons.AddOption(	'--mkl',
				dest='mkl',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Enable Intel MKL [default: %default]"
		)
		self.mkl = scons.GetOption('mkl')


		scons.AddOption(	'--mic',
				dest='mic',
				type='choice',
				choices=['enable', 'disable'],
				default='disable',
				help="Enable Intel MIC (XeonPhi) [default: %default]"
		)
		self.mic = scons.GetOption('mic')


		#
		# LIB XML
		#
		scons.AddOption(	'--libxml',
				dest='libxml',
				type='choice',
				choices=['enable','disable'],
				default='disable',
				help='Compile with libXML related code: enable, disable [default: %default]'
		)
		self.libxml = scons.GetOption('libxml')


		scons.AddOption(	'--plane-spectral-dealiasing',
				dest='plane_spectral_dealiasing',
				type='choice',
				choices=['enable','disable'],
				default='enable',
				help='spectral dealiasing (3N/2-1 rule): enable, disable [default: %default]'
		)
		self.plane_spectral_dealiasing = scons.GetOption('plane_spectral_dealiasing')

		scons.AddOption(	'--sphere-spectral-dealiasing',
				dest='sphere_spectral_dealiasing',
				type='choice',
				choices=['enable','disable'],
				default='enable',
				help='spectral dealiasing (3N/2-1 rule): enable, disable [default: %default]'
		)
		self.sphere_spectral_dealiasing = scons.GetOption('sphere_spectral_dealiasing')

		# Always use dealiasing for sphere
		if self.sphere_spectral_space == 'enable' and self.sphere_spectral_dealiasing != 'enable':
			raise Exception("self.sphere_spectral_dealiasing != enable")

		scons.AddOption(	'--quadmath',
				dest='quadmath',
				type='choice',
				choices=['enable','disable'],
				default='disable',
				help='quadmath: enable, disable [default: %default]'
		)
		self.quadmath = scons.GetOption('quadmath')



		scons.AddOption(	'--gui',
				dest='gui',
				type='choice',
				choices=['enable','disable'],
				default='disable',
				help='gui: enable, disable [default: %default]'
		)
		self.gui = scons.GetOption('gui')




		scons.AddOption(	'--rexi-thread-parallel-sum',
				dest='rexi_thread_parallel_sum',
				type='choice',
				choices=['enable','disable'],
				default='disable',
				help='Use a par for loop over the sum in REXI: (enable, disable) [default: %default]\n\tWARNING: This also disables the parallelization-in-space with OpenMP'
		)
		self.rexi_thread_parallel_sum = scons.GetOption('rexi_thread_parallel_sum')


		scons.AddOption(	'--rexi-timings',
				dest='rexi_timings',
				type='choice',
				choices=['enable','disable'],
				default='disable',
				help='REXI timings: enable, disable [default: %default]'
		)
		self.rexi_timings = scons.GetOption('rexi_timings')

		scons.AddOption(	'--rexi-timings-additional-barriers',
				dest='rexi_timings_additional_barriers',
				type='choice',
				choices=['enable','disable'],
				default='disable',
				help='REXI timings with additional barriers: enable, disable [default: %default]\nThis is helpful for improved measurements if TurboBoost is activated'
		)
		self.rexi_timings_additional_barriers = scons.GetOption('rexi_timings_additional_barriers')


		scons.AddOption(	'--sweet-mpi',
				dest='sweet_mpi',
				type='choice',
				choices=['enable','disable'],
				default='disable',
				help='Enable MPI commands e.g. for REXI sum: enable, disable [default: %default]'
		)
		self.sweet_mpi = scons.GetOption('sweet_mpi')



		scons.AddOption(	'--parareal',
				dest='parareal',
				type='choice',
				choices=['none', 'serial','mpi'],
				default='none',
				help='Enable Parareal (none, serial, mpi) [default: %default]\nOnly works, if Parareal is supported by the simulation'
		)
		self.parareal = scons.GetOption('parareal')


		files = os.listdir('src/programs/')
		files = sorted(files)
		example_programs = []
		for f in files:
			if os.path.isfile('src/programs/'+f):
				example_programs.append(f[0:-4])

		scons.AddOption(      '--program',
				dest='program',
				type='choice',
				choices=example_programs,
				default='',
				help='Specify program to compile: '+', '.join(example_programs)+' '*80+' [default: %default]'
		)
		self.program = scons.GetOption('program')


		scons.AddOption(      '--program-binary-name',
				dest='program_binary_name',
				type='string',
				action='store',
				help='Name of program binary, default: [auto]',
				default=''
		)
		self.program_binary_name = scons.GetOption('program_binary_name')



		files = os.listdir('src/unit_tests/')
		files = sorted(files)
		unit_tests_programs = []
		for f in files:
			unit_tests_programs.append(f[0:-4])


		scons.AddOption(      '--unit-test',
				dest='unit_test',
				type='choice',
				choices=unit_tests_programs,
				default='',
				help='Specify unit tests to compile: '+', '.join(unit_tests_programs)+' '*80+' [default: %default]'
		)
		self.unit_test = scons.GetOption('unit_test')


		threading_constraints = ['off', 'omp']
		scons.AddOption(	'--threading',
				dest='threading',
				type='choice',
				choices=threading_constraints,
				default='omp',
				help='Threading to use '+' / '.join(threading_constraints)+', default: off'
		)
		self.threading = scons.GetOption('threading')


		scons.AddOption(	'--cxx-flags',
				dest='cxx_flags',
				type='string',	
				default='',
				help='Additional cxx-flags, default: ""'
		)

		self.cxx_flags = scons.GetOption('cxx_flags')


		scons.AddOption(	'--ld-flags',
				dest='ld_flags',
				type='string',	
				default='',
				help='Additional ld-flags, default: ""'
		)



	def getProgramName(self, ignore_errors = False):
		self.makeOptionsConsistent()

		if self.program != '':
			self.program_name = self.program

		elif self.unit_test != '':
			self.program_name = self.unit_test

		else:
			self.program_name = 'DUMMY'
			self.info("")
			self.info("")
			self.info("Neither a program name, nor a unit test is given:")
			self.info("  use --program=[program name] to specify the program")
			self.info("  or --unit-test=[unit test] to specify a unit test")
			self.info("")
			if not ignore_errors:
				sys.exit(1)


		if self.program_binary_name != '':
			return self.program_binary_name

		exec_name = self.program_name

		if self.plane_spectral_space == 'enable':
			exec_name+='_plspec'

		if self.plane_spectral_dealiasing == 'enable':
			exec_name+='_pldeal'

		if self.sphere_spectral_space == 'enable':
			exec_name+='_spspec'

		if self.sphere_spectral_dealiasing == 'enable':
			exec_name+='_spdeal'

		if self.gui == 'enable':
			exec_name+='_gui'

		if self.quadmath == 'enable':
			exec_name+='_quadmath'

		if self.sweet_mpi == 'enable':
			exec_name+='_mpi'

		if self.threading in ['omp']:
			exec_name+='_'+self.threading
		else:
			if self.llvm_omp_override:
				print("WARNING: adding _omp despite program was not compiled with OpenMP activated. This is for compatibility reasons only!")
				exec_name+='_omp'
			
		if self.rexi_thread_parallel_sum == 'enable':
			exec_name+='_rxthpar'

		if self.rexi_timings == 'enable':
			exec_name+='_rxtime'

		if self.rexi_timings_additional_barriers == 'enable':
			exec_name+='_rxtbar'


		if self.numa_block_allocator in [1, 2]:
			exec_name+='_numa'+str(self.numa_block_allocator)

		if self.libfft == 'enable':
			exec_name+='_fft'

		if self.llvm_gnu_override:
			print("WARNING: adding _omp despite program was not compiled with LLVM. This is for compatibility reasons only!")
			exec_name += '_gnu'
		else:
			exec_name += '_'+self.compiler

		exec_name += '_'+self.mode

		return exec_name


	def getUniqueID(self, i_filter_list = []):
		"""
		Return a unique ID including *all* string and number attributes of this class

		i_filter:
			list of filter names to filter out from unique ID generation
		"""
		return self.getProgramName()



	def getUniqueParID(self, i_filter_list = []):
		"""
		TODO: Get rid of this!
		"""
		if not 'compiler_parallelization' in i_filter_list:
			retval = 'COMP'
			if self.sweet_mpi == 'enable':
				retval+='_mpi'

			if self.threading in ['omp']:
				retval+='_'+self.threading

			if self.rexi_thread_parallel_sum == 'enable':
				retval+='_rxthpar'

		return retval



	def getOptionList(self):
		return self.__dict__


if __name__ == "__main__":

	p = SWEETCompileOptions()

	s = p.getSConsParams()
	p.info(s)

	s = p.getProgramName(True)
	p.info(s)

	p.info("FIN")
