import os
import subprocess
import re
import sys
import platform


sys.path.append("./python_mods")
from SWEETCompileOptions import *
sys.path.remove("./python_mods")

#
# Setup parallel compilation
#
import multiprocessing
num_cpu = multiprocessing.cpu_count()
SetOption('num_jobs', num_cpu)


def exec_command(command):
	process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	# combine stdout and stderr
	out = out+err
	out = out.decode("utf-8")
	out = out.replace("\r", "")
	return out

#
# determine hostname
#
hostname = exec_command('hostname')
hostname = hostname.replace("\n", "")

env = Environment(ENV = os.environ)

if 'SWEET_ROOT' not in os.environ:
	print("*"*80)
	print("* Welcome to SWEET, Awesome!")
	print("*"*80)
	print("* The SWEET_ROOT environment variable was not found.")
	print("* Please load all SWEET environment variables with")
	print("*   $ source ./local_software/env_vars.sh")
	print("* or")
	print("*   $ . ./local_software/env_vars.sh")
	print("* (including the single dot at the beginning)")
	print("*"*80)
	Exit(1)

env['SWEET_ROOT'] = os.environ['SWEET_ROOT']

p = SWEETCompileOptions()



###################################################################
# fix LD LIB PATH
###################################################################
if 'LD_LIBRARY_PATH' in os.environ:
	env.Append(LIBPATH=os.environ['LD_LIBRARY_PATH'].split(':'))



###########################################
# Compile options
###########################################


p.sconsProcessOptions()

p.makeOptionsConsistent()

if p.libfft == 'enable':

	env.Append(CXXFLAGS = ' -DSWEET_USE_LIBFFT=1')

	if p.mkl == 'enable':
		print("INFO: Using Intel MKL instead of FFTW");

		if p.threading != 'omp':
			env.Append(CXXFLAGS=['-mkl=sequential'])
			env.Append(LINKFLAGS=['-mkl=sequential'])
		else:
			env.Append(CXXFLAGS=['-mkl=parallel'])
			env.Append(LINKFLAGS=['-mkl=parallel'])

	else:
		if p.threading == 'omp':
			env.Append(LIBS=['fftw3_omp'])

		env.Append(LIBS=['fftw3'])
else:
	env.Append(CXXFLAGS = ' -DSWEET_USE_LIBFFT=0')



if p.libxml == 'enable':
	env.ParseConfig("xml2-config --cflags --libs")

if p.parareal == 'mpi':
	raise Exception("TODO: Implement MPI Parareal")

p.ld_flags = GetOption('ld_flags')

env.Append(CXXFLAGS=' -DSWEET_SIMD_ENABLE='+('1' if p.simd=='enable' else '0'))
env.Append(CXXFLAGS=' -DCONFIG_ENABLE_LIBXML='+('1' if p.libxml=='enable' else '0'))
env.Append(CXXFLAGS = p.cxx_flags)
env.Append(LINKFLAGS = p.ld_flags)



p.llvm_gnu_override = False
p.llvm_omp_override = False


if p.plane_spectral_space == 'enable':
	env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_SPACE=1')
else:
	env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_SPACE=0')



if p.plane_spectral_dealiasing == 'enable':
	env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_DEALIASING=1')
else:
	env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_DEALIASING=0')




if p.sphere_spectral_space == 'enable':

#	if p.numa_block_allocator != 0:
#		print("NUMA Block allocator option > 0 not tested with SHTNS, buffers in SHTNS may be shared")
#		Exit(1)

	env.Append(CXXFLAGS = ' -DSWEET_USE_SPHERE_SPECTRAL_SPACE=1')
else:
	env.Append(CXXFLAGS = ' -DSWEET_USE_SPHERE_SPECTRAL_SPACE=0')


if p.sphere_spectral_dealiasing == 'enable':
	env.Append(CXXFLAGS = ' -DSWEET_USE_SPHERE_SPECTRAL_DEALIASING=1')
else:
	if p.sphere_spectral_space == 'enable':
		raise Exception("No anti-aliasing on sphere as compile option supported, please simply use command line options to specify lower physical resolution!")

	env.Append(CXXFLAGS = ' -DSWEET_USE_SPHERE_SPECTRAL_DEALIASING=0')






if p.parareal == 'none':
	env.Append(CXXFLAGS = ' -DSWEET_PARAREAL=0')
elif p.parareal == 'serial':
	env.Append(CXXFLAGS = ' -DSWEET_PARAREAL=1')
elif p.parareal == 'mpi':
	env.Append(CXXFLAGS = ' -DSWEET_PARAREAL=2')
else:
	print("Invalid option '"+str(p.parareal)+"' for parareal method")
	sys.exit(1)
	



if p.compiler == 'gnu':
	reqversion = [4,6,1]

	#
	# get gcc version using -v instead of -dumpversion since SUSE gnu compiler
	# returns only 2 instead of 3 digits with -dumpversion
	#
	gccv = exec_command('g++ -v').splitlines()

	# updated to search for 'gcc version ' line prefix
	found = False
	found_line = ''

	search_string = 'gcc version '
	for l in gccv:
		if l[:len(search_string)] == search_string:
			found_line = l
			found = True

			gccversion = found_line.split(' ')[2].split('.')
			break

	if not found:
		print(search_string+" not found... testing for next search string")
		search_string = 'gcc-Version '
		for l in gccv:
			if l[:len(search_string)] == search_string:
				found_line = l
				found = True
	
				gccversion = found_line.split(' ')[1].split('.')
				break

	if not found:
		print(search_string+" not found... testing if this is LLVM on MacOSX")
		found = False
		for l in gccv:
			if 'Apple LLVM' in l:
				found = True
				break
		if not found:
			print("LLVM not detected")
			sys.exit(1)

		p.llvm_gnu_override = True
		p.compiler = 'llvm'

	else:
		for i in range(0, 3):
			if (int(gccversion[i]) > int(reqversion[i])):
				break
			if (int(gccversion[i]) < int(reqversion[i])):
				print('At least GCC Version 4.6.1 necessary.')
				Exit(1)

	if p.compiler == 'gnu':
		# eclipse specific flag
		env.Append(CXXFLAGS=' -fmessage-length=0')

		# c++0x flag
		env.Append(CXXFLAGS=' -std=c++0x')

		# be pedantic to avoid stupid programming errors
	#	env.Append(CXXFLAGS=' -pedantic')

		# speedup compilation - remove this when compiler slows down or segfaults by running out of memory
		env.Append(CXXFLAGS=' -pipe')

		# activate gnu C++ compiler

		if p.fortran_source=='enable':
			#env.Replace(F90='gfortran')
			env.Append(F90FLAGS=' -cpp')
			env.Append(LIBS=['gfortran'])

	#	env.Replace(CXX = 'g++-4.7')
	#	env.Replace(CXX = 'g++')



if p.compiler == 'intel':
	reqversion = [12,1]
	iccversion_line = exec_command('icpc -dumpversion -w').splitlines()[0]

	if iccversion_line != 'Mainline':
		iccversion = iccversion_line.split('.')
		for i in range(0, 2):
			if (int(iccversion[i]) > int(reqversion[i])):
				break
			if (int(iccversion[i]) < int(reqversion[i])):
				print('ICPC Version 12.1 necessary.')
				Exit(1)

	# override g++ called by intel compiler to determine location of header files
	if p.gxx_toolchain != '':
		env.Append(CXXFLAGS='-gxx-name='+p.gxx_toolchain)
		env.Append(LINKFLAGS='-gxx-name='+p.gxx_toolchain)

	env.Append(LINKFLAGS='-shared-intel')
	env.Append(LINKFLAGS='-shared-libgcc')
	env.Append(LINKFLAGS='-debug inline-debug-info')

	# eclipse specific flag
	env.Append(CXXFLAGS=' -fmessage-length=0')

	# c++0x flag
	env.Append(CXXFLAGS=' -std=c++0x')

	# output more warnings
	env.Append(CXXFLAGS=' -w1')


	# compiler option which has to be appended for icpc 12.1 without update1
#	env.Append(CXXFLAGS=' -U__GXX_EXPERIMENTAL_CXX0X__')


	# UBUNTU FIX for i386 systems
	lines = exec_command('uname -i').splitlines()

	for i in lines:
		if i == 'i386':
			env.Append(CXXFLAGS=' -I/usr/include/i386-linux-gnu/')

	# SSE 4.2
	#env.Replace(CXX = 'icpc')
	#env.Replace(LINK='icpc')

	if p.fortran_source == 'enable':
		env.Append(LIBS=['gfortran'])
		env.Append(LIBS=['ifcore'])
		#env.Replace(F90='ifort')
		env.Append(F90FLAGS=' -fpp')



# WARNING: don't use 'elif' here wince llvm may be activated via the 'gnu' compiler option
if p.compiler == 'llvm':
	reqversion = [3,1]

	if p.gxx_toolchain != '':
		env.Append(CXXFLAGS=' --gcc-toolchain='+p.gxx_toolchain)

	# eclipse specific flag
	env.Append(CXXFLAGS=' -fmessage-length=0')

	# c++0x flag
	env.Append(CXXFLAGS=' -std=c++0x')

	# support __float128
#	env.Append(CXXFLAGS=' -D__STRICT_ANSI__')

	# be pedantic to avoid stupid programming errors
	env.Append(CXXFLAGS=' -pedantic')

	# speedup compilation - remove this when compiler slows down or segfaults by running out of memory
	env.Append(CXXFLAGS=' -pipe')

	# activate gnu C++ compiler


	if p.sphere_spectral_space == 'enable':
		# append gfortran library
		env.Append(LIBS=['gfortran'])

	if p.threading == 'omp':
		print("")
		print('WARNING: OpenMP with LLVM not supported, deactivating')
		print("")

		p.llvm_omp_override = True
		p.threading = 'off'

	# todo: fix me also for intel mpicxx compiler
	#env.Replace(CXX = 'clang++')

	if p.fortran_source == 'enable':
		env.Append(LIBS=['gfortran'])
		print("TODO: LLVM compiler not yet supported with fortran enabled")
		Exit(-1)



if p.mode == 'debug':
	env.Append(CXXFLAGS=' -DSWEET_DEBUG=1')

	if p.compiler == 'gnu':
		env.Append(CXXFLAGS='-O0 -g3 -Wall')

		# integer overflow check
		env.Append(CXXFLAGS=' -ftrapv')

	elif p.compiler == 'llvm':
		env.Append(CXXFLAGS='-O0 -g3 -Wall')

	elif p.compiler == 'intel':
		env.Append(CXXFLAGS='-O0 -g -traceback')
#		env.Append(CXXFLAGS=' -fp-trap=common')

	elif p.compiler == 'pgi':
		env.Append(CXXFLAGS='-O0 -g -traceback')


	if p.fortran_source == 'enable':
		if p.compiler == 'gnu':
			env.Append(F90FLAGS='-g -O0')
		elif p.compiler == 'intel':
			env.Append(F90FLAGS='-g -O0 -traceback')


elif p.mode == 'release':
	env.Append(CXXFLAGS=' -DSWEET_DEBUG=0')

	# deactivate assertion calls
	env.Append(CXXFLAGS=' -DNDEBUG=1')

	if p.compiler == 'gnu':
		env.Append(CXXFLAGS=' -O3 -mtune=native')

	elif p.compiler == 'llvm':
		env.Append(CXXFLAGS=' -O3 -mtune=native')

	elif p.compiler == 'intel':
		env.Append(CXXFLAGS=' -O2 -fno-alias')

		if p.mic != 'enable':
			env.Append(CXXFLAGS=' -xHost')

	elif p.compiler == 'pgi':
		env.Append(CXXFLAGS='-O3 -fast -Mipa=fast,inline -Msmartalloc')

	if p.fortran_source == 'enable':
		if p.compiler == 'gnu':
			env.Append(F90FLAGS=' -O2')
		elif p.compiler == 'intel':
			env.Append(F90FLAGS=' -O2')


if p.quadmath == 'enable':
	env.Append(CXXFLAGS=' -DSWEET_QUADMATH=1')
	env.Append(LIBS=['quadmath'])
else:
	env.Append(CXXFLAGS=' -DSWEET_QUADMATH=0')



if p.gui == 'enable':
	# compile flags
	env.Append(CXXFLAGS=' -I'+os.environ['HOME']+'/local/include')
	env.Append(CXXFLAGS=' -DSWEET_GUI=1')


	# linker flags

	if exec_command('uname -s') == "Darwin":
		# ASSUME MACOSX SYSTEM
		env.Append(LINKFLAGS='-framework OpenGL')
	else:
		env.Append(LIBS=['GL'])

	reqversion = [2,0,0]
	sdlversion = exec_command('sdl2-config --version')
	sdlversion = sdlversion.replace("\n", "").split('.')

	for i in range(0, 3):
		if (int(sdlversion[i]) > int(reqversion[i])):
			break;
		if (int(sdlversion[i]) < int(reqversion[i])):
			print('libSDL Version 2.0.0 necessary.')
			Exit(1)

	env.ParseConfig("sdl2-config --cflags --libs")
	env.ParseConfig("freetype-config --cflags --libs")
else:
	env.Append(CXXFLAGS=' -DSWEET_GUI=0')

  

if p.fortran_source == 'enable':
	env.Append(CXXFLAGS=' -DSWEET_FORTRAN=1')
else:
	env.Append(CXXFLAGS=' -DSWEET_FORTRAN=0')


env.Append(CXXFLAGS=' -DNUMA_BLOCK_ALLOCATOR_TYPE='+str(p.numa_block_allocator))

if p.numa_block_allocator in [1, 2]:
	env.Append(LIBS=['numa'])


if p.threading == 'omp':
	env.Append(CXXFLAGS=['-fopenmp'])
	env.Append(LINKFLAGS=['-fopenmp'])

	env.Append(CXXFLAGS=' -DSWEET_SPACE_THREADING=1')
else:
	env.Append(CXXFLAGS=' -DSWEET_SPACE_THREADING=0')


if p.pfasst_cpp == 'enable':
	env.Append(CXXFLAGS=['-Ilocal_software/local/include/eigen3'])
	env.Append(CXXFLAGS=['-DSWEET_PFASST_CPP=1'])
else:
	env.Append(CXXFLAGS=['-DSWEET_PFASST_CPP=0'])

if p.eigen == 'enable':
	env.Append(CXXFLAGS=['-Ilocal_software/local/include/eigen3'])
	env.Append(CXXFLAGS=['-DSWEET_EIGEN=1'])
else:
	env.Append(CXXFLAGS=['-DSWEET_EIGEN=0'])


if p.libpfasst == 'enable':
	env.Append(CXXFLAGS=['-Llibpfasst'])
	env.Append(CXXFLAGS=['-DSWEET_LIBPFASST=1'])

	env.Append(LIBS=['libpfasst'])

	# enable MPI per default for libpfasst
	p.sweet_mpi = 'enable'

else:
	env.Append(CXXFLAGS=['-DSWEET_LIBPFASST=0'])




if p.sweet_mpi == 'enable':
	env.Append(CXXFLAGS = ' -DSWEET_MPI=1')
else:
	env.Append(CXXFLAGS = ' -DSWEET_MPI=0')


if p.sweet_mpi == 'enable':
	print("Enabling MPI for REXI")
	print("Warning: Compiler checks not done")

	if p.compiler == 'gnu':
		#env.Replace(CXX = 'mpiCC')
		#env.Replace(LINK = 'mpiCC')
		#env.Replace(F90 = 'mpif90')

		# GNU compiler needs special treatment!
		# Linking with Fortran MPI requires
		# for OpenMPI: -lmpi_mpifh
		# for MPICH: -lmpif90

		output = exec_command('mpiCC -v')
		if 'MPICH' in output:
			if p.fortran_source == 'enabled':
				env.Append(LINKFLAGS='-lmpif90')
		else:
			pass
			#print("*"*80)
			#print("ERROR: MPI VERSION NOT DETECTED, ASSUMING OPENMPI!!!!!")
			#print("*"*80)
			#env.Append(LINKFLAGS='-lmpi_mpifh')
			
		

	elif p.compiler == 'intel':
		#env.Replace(CXX = 'mpiicpc')
		#env.Replace(LINK = 'mpiicpc')
		#env.Replace(F90 = 'mpif90')
		pass

	if p.threading != 'off' and p.compiler == 'intel':
		env.Append(CXXFLAGS='-mt_mpi')
		env.Append(LINKFLAGS='-mt_mpi')



#
# Override compiler settings from environment variable
#
override_list = ['CC', 'CXX', 'F90', 'LINK']
for i in override_list:
	if p.sweet_mpi == 'enable':
		if 'SWEET_MPI'+i in env['ENV']:
			print("INFO: Overriding environment variable "+i+"="+env['ENV']['SWEET_MPI'+i])
			env[i] = env['ENV']['SWEET_MPI'+i]


	else:
		if 'SWEET_'+i in env['ENV']:
			print("INFO: Overriding environment variable "+i+"="+env['ENV']['SWEET_'+i])
			env[i] = env['ENV']['SWEET_'+i]



if p.libsph == 'enable':
	if p.threading == 'omp':
		env.Append(LIBS=['shtns_omp'])
	else:
		env.Append(LIBS=['shtns'])

	# Add LAPACK libraries
	if p.lapack == 'enable':
		env.Append(LIBS=['lapack'])
		env.Append(LIBS=['blas'])
		env.Append(CXXFLAGS=['-DSWEET_LAPACK=1'])
	else:
		env.Append(CXXFLAGS=['-DSWEET_LAPACK=0'])

	if p.compiler == 'gnu':
		env.Append(LIBS=['gfortran'])



if p.mic == 'enable':
	env.Append(CXXFLAGS=['-mmic'])
	env.Append(LINKFLAGS=['-mmic'])


if p.rexi_thread_parallel_sum == 'enable' and p.threading == 'omp':
	raise Exception('ERROR: "REXI Parallel Sum" and "Threading" is both activated')



#
# If SWEET_REXI_THREAD_PARALLEL_SUM is activated, the REXI sum is computed
# with parallel for over the sum terms
#
if p.rexi_thread_parallel_sum == 'enable':
	# Same for gcc/icpc
	env.Append(LINKFLAGS=['-fopenmp'])

	#
	# Also add CXX flags here.
	# This requires that *ALL* space-related 
	#
	#pragma omp parallel for
	#
	# loops are annoted with
	#
	#if SWEET_SPACE_THREADING
	#
	env.Append(CXXFLAGS=['-fopenmp'])

	# Compile flag is set in sconscript

	# Activate precompiler flag
	env.Append(CXXFLAGS=' -DSWEET_REXI_THREAD_PARALLEL_SUM=1')
else:
	env.Append(CXXFLAGS=' -DSWEET_REXI_THREAD_PARALLEL_SUM=0')


if p.rexi_timings == 'enable':
	env.Append(CXXFLAGS=' -DSWEET_REXI_TIMINGS=1')
else:
	env.Append(CXXFLAGS=' -DSWEET_REXI_TIMINGS=0')

if p.rexi_timings_additional_barriers == 'enable':
	env.Append(CXXFLAGS=' -DSWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS=1')
else:
	env.Append(CXXFLAGS=' -DSWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS=0')


if p.debug_symbols == 'enable':
	env.Append(CXXFLAGS = '-g')
	env.Append(LINKFLAGS = '-g')

	if p.compiler == 'intel':
		env.Append(CXXFLAGS = ' -shared-intel -shared-libgcc -debug inline-debug-info')
		env.Append(LINKFLAGS = ' -shared-intel  -shared-libgcc -debug inline-debug-info')


exec_name = p.getProgramName()


#
# BUILD directory
#
build_dir='/tmp/scons_build_'+exec_name+'/'

if p.libpfasst == 'enable':
	#env.Append(F90FLAGS = ['-Ilocal_software/local_src/libpfasst/include'])
	env.Append(F90FLAGS = ['-Ilocal_software/local/include'])
#
# USE build directory for Fortran module output
#
env.Append(F90FLAGS = '-J'+build_dir)


env.Append(CPPPATH = ['/usr/local/include', '/usr/include'])




#
# FORTRAN stuff
#

#fortran_mod_dir = build_dir+'/fortran_mods'



#if not os.path.exists(fortran_mod_dir):
#    os.makedirs(fortran_mod_dir)
#env = env.Clone(FORTRANMODDIR = fortran_mod_dir)


# also include the 'src' directory to search for dependencies
env.Append(CPPPATH = ['.', './src/', './src/include'])
# also for Fortran!
env.Append(F90PATH = ['.', './src/'])



######################
# get source code files
######################

env.src_files = []

# local software directories
env.Append(LINKFLAGS=['-L./local_software/local/lib'])
env.Append(CPPPATH=['./local_software/local/include'])

if p.program_name != 'DUMMY':

	env.SConscript('sconscript', variant_dir=build_dir, duplicate=0, exports=['env', 'p'])

	print('')
	print('            Program: '+p.program_name)
	print('Building executable: '+exec_name)
	print('')

	obj_files = []
	for i in env.src_files:
		for a in i:
			if str(a).split(".")[-1] == "o":
				obj_files.append(str(a))

	env.Program('build/'+exec_name, obj_files)

