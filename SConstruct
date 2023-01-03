import os
import subprocess
import re
import sys
import platform


print(" ".join(sys.argv))


# SCons doesn't inherit the python environment variables
sys.path.append("./mule_local/python/mule_local")
from JobCompileOptions import *
sys.path.remove("./mule_local/python/mule_local")

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
# uname information
#
uname = exec_command('uname').replace("\n", "")


#
# determine hostname
#
hostname = exec_command('hostname')
hostname = hostname.replace("\n", "")

env = Environment(ENV = os.environ)

if 'MULE_SOFTWARE_ROOT' not in os.environ:
    print("*"*80)
    print("* Welcome to SWEET, Awesome!")
    print("*"*80)
    print("* The MULE_SOFTWARE_ROOT environment variable was not found.")
    print("* Please load all SWEET environment variables with")
    print("*   $ source ./local_software/env_vars.sh")
    print("* or")
    print("*   $ . ./local_software/env_vars.sh")
    print("* (including the single dot at the beginning)")
    print("*"*80)
    Exit(1)

env['MULE_SOFTWARE_ROOT'] = os.environ['MULE_SOFTWARE_ROOT']

p = JobCompileOptions()


###################################################################
# fix LD LIB PATH
###################################################################
if 'LD_LIBRARY_PATH' in os.environ:
    env.Append(LIBPATH=os.environ['LD_LIBRARY_PATH'].split(':'))



###########################################
# Compile options
###########################################


p.sconsProcessOptions()

#
# Process scons options provided within source file
# e.g. 
# * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere/
# * MULE_COMPILE_FILES_AND_DIRS: src/include/benchmarks_sphere_swe/
# * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
#
pso_ = p.get_program_specific_options()
p.process_scons_options(pso_['scons_options'])

#
# Cleanup options
#
p.makeOptionsConsistent()



if p.libxml == 'enable':
    env.ParseConfig("xml2-config --cflags --libs")

####if p.parareal == 'mpi':
####    raise Exception("TODO: Implement MPI Parareal")

p.ld_flags = GetOption('ld_flags')

env.Append(CXXFLAGS=' -DSWEET_SIMD_ENABLE='+('1' if p.simd=='enable' else '0'))
env.Append(CXXFLAGS=' -DCONFIG_ENABLE_LIBXML='+('1' if p.libxml=='enable' else '0'))
env.Append(CXXFLAGS = p.cxx_flags)
env.Append(LINKFLAGS = p.ld_flags)




if p.plane_spectral_space == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_SPACE=1')
else:
    env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_SPACE=0')



if p.plane_spectral_dealiasing == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_DEALIASING=1')
else:
    env.Append(CXXFLAGS = ' -DSWEET_USE_PLANE_SPECTRAL_DEALIASING=0')




if p.sphere_spectral_space == 'enable':
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
    raise Exception("Invalid option '"+str(p.parareal)+"' for parareal method")


if p.parareal_scalar == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_PARAREAL_SCALAR=1')
if p.parareal_plane == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_PARAREAL_PLANE=1')
if p.parareal_sphere == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_PARAREAL_SPHERE=1')
if p.parareal_plane_swe == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_PARAREAL_PLANE_SWE=1')
if p.parareal_plane_burgers == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_PARAREAL_PLANE_BURGERS=1')


if p.xbraid == 'none':
    env.Append(CXXFLAGS = ' -DSWEET_XBRAID=0')
elif p.xbraid == 'mpi':
    env.Append(CXXFLAGS=['-Ilocal_software/local/include/xbraid'])
    env.Append(LIBS=['braid'])
    env.Append(CXXFLAGS = ' -DSWEET_XBRAID=1')
else:
    raise Exception("Invalid option '"+str(p.xbraid)+"' for XBraid")


if p.xbraid_scalar == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_XBRAID_SCALAR=1')
if p.xbraid_plane == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_XBRAID_PLANE=1')
if p.xbraid_sphere == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_XBRAID_SPHERE=1')
if p.xbraid_plane_swe == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_XBRAID_PLANE_SWE=1')
if p.xbraid_plane_burgers == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_XBRAID_PLANE_BURGERS=1')

if p.scalar_type == 'double':
    env.Append(CXXFLAGS = ' -DSWEET_SCALAR_COMPLEX=0')
elif p.scalar_type == 'complex':
    env.Append(CXXFLAGS = ' -DSWEET_SCALAR_COMPLEX=1')
else:
    raise Exception("Invalid option '"+str(p.scalar_type)+"' for scalar_type")

if int(p.N_ode) > 0:
    env.Append(CXXFLAGS = ' -DSWEET_N_ODE='+str(p.N_ode))

#
# Override compiler settings from environment variable
#
override_list = [
        'CC', 'CCFLAGS',
        'CXX', 'CXXFLAGS',
        'F90', 'F90FLAGS',
        'FC', 'FCFLAGS',
        'LINK', 'LINKFLAGS',
        'LIBS',
        #'LD'
        ]
for i in override_list:

    if p.sweet_mpi == 'enable':
        mi = 'MULE_MPI'+i
        if mi in env['ENV']:
            if 'FLAGS' in mi:
                print("INFO: Appending to "+i+"+= "+env['ENV'][mi])
                env[i].append(env['ENV'][mi])
            elif 'LIBS'==i:
                print("INFO: Appending to "+i+"+= "+env['ENV'][mi])
                env.Append(LIBS=env['ENV'][mi])
            else:
                print("INFO: Using MULE_MPI* environment variable to set "+i+"="+env['ENV'][mi])
                env[i] = env['ENV'][mi]
            continue

    if i in env['ENV']:
        if 'FLAGS' in i:
            print("INFO: Appending to "+i+"+= "+env['ENV'][i])
            env[i] = env[i]+' '+env['ENV'][i]
        if 'LIBS' == i:
            print("INFO: Appending to "+i+"+= "+env['ENV'][i])
            env.Append(LIBS=env['ENV'][i])
        else:
            print("INFO: Overriding environment variable "+i+"="+env['ENV'][i])
            env[i] = env['ENV'][i]



#
# This variable is used to select the compiler which is used.
# Note, that this might not match the compiler specified in the job ID!
#
compiler_to_use = p.compiler

#
# Autodetect compiler for MPI
#
if p.sweet_mpi == 'enable':
    env.Append(CXXFLAGS = ' -DSWEET_MPI=1')
    mpicc_version = exec_command(os.environ['MULE_MPICC']+' --version')

    if mpicc_version[0:3] == 'icc' or mpicc_version[0:4] == 'icpc':
        if p.compiler != 'intel':
            print("Autodetected Intel compiler for MPI, switching to Intel compiler options")
            compiler_to_use = 'intel'

    elif mpicc_version[0:3] == 'gcc' or mpicc_version[0:3] == 'g++':
        if p.compiler != 'gnu':
            print("Autodetected GNU compiler for MPI, switching to GNUl compiler options")
            compiler_to_use = 'gnu'

else:
    env.Append(CXXFLAGS = ' -DSWEET_MPI=0')


if compiler_to_use == 'gnu':

    if False:
        reqversion = [5,0,0]

        #
        # get gcc version using -v instead of -dumpversion since SUSE gnu compiler
        # returns only 2 instead of 3 digits with -dumpversion
        #

        gccv = exec_command(env['CC']+' -v').splitlines()

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

            p.compiler = 'llvm'

        else:
            for i in range(0, 3):
                if (int(gccversion[i]) > int(reqversion[i])):
                    break
                if (int(gccversion[i]) < int(reqversion[i])):
                    print('At least GCC Version 4.6.1 necessary.')
                    Exit(1)

    # eclipse specific flag
    env.Append(CXXFLAGS=' -fmessage-length=0')

    # c++0x flag
    env.Append(CXXFLAGS=' -std=c++0x')

    # be pedantic to avoid stupid programming errors
    #  env.Append(CXXFLAGS=' -pedantic')

    # speedup compilation - remove this when compiler slows down or segfaults by running out of memory
    #env.Append(CXXFLAGS=' -pipe')

    if p.fortran_source == 'enable':
        #env.Replace(F90='gfortran')
        env.Append(F90FLAGS=' -cpp')
        env.Append(LIBS=['gfortran'])



if p.libpfasst == 'enable':
    env.Append(CXXFLAGS=['-Llibpfasst'])
    env.Append(CXXFLAGS=['-DSWEET_LIBPFASST=1'])

    env.Append(LIBS=['libpfasst'])

    # enable MPI per default for libpfasst
    p.sweet_mpi = 'enable'

else:
    env.Append(CXXFLAGS=['-DSWEET_LIBPFASST=0'])



env.Append(LIBS=['m'])

if p.sweet_mpi == 'enable':
    print("Enabling MPI for REXI")
    print("Warning: Compiler checks not done")

    if compiler_to_use == 'gnu':
        # GNU compiler needs special treatment!
        # Linking with Fortran MPI requires
        # for OpenMPI: -lmpi_mpifh
        # for MPICH: -lmpifort

        output = exec_command(env['CC']+' -v')
        if 'MPICH' in output:
            if p.fortran_source == 'enable':
                env.Append(LIBS='mpifort')
        else:
            # Assume OpenMPI
            if p.fortran_source == 'enable':
                env.Append(LIBS='mpi_cxx')
        
    if p.threading != 'off' and compiler_to_use == 'intel':
        env.Append(CXXFLAGS='-mt_mpi')
        env.Append(LINKFLAGS='-mt_mpi')



if compiler_to_use == 'intel':
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
#    env.Append(CXXFLAGS=' -U__GXX_EXPERIMENTAL_CXX0X__')


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



# WARNING: don't use 'elif' here since clang may be activated via the 'gnu' compiler option
if compiler_to_use == 'llvm':
    reqversion = [9,0]
    if p.threading == 'omp':
        reqversion = [9,0]
    version_line = exec_command('clang++ --version').splitlines()[0]

    verpos = version_line.find(" version ")
    if verpos == -1:
        raise Exception("Version number not detected")

    version_line = version_line[verpos+9:]
    version = version_line.split('.')

    for i in range(0, 2):
        if (int(version[i]) > int(reqversion[i])):
            break
        if (int(version[i]) < int(reqversion[i])):
            print("LLVM version "+(".".join(reqversion))+" or higher required")
            Exit(1)

    #if p.gxx_toolchain != '':
    #    env.Append(CXXFLAGS=' --gcc-toolchain='+p.gxx_toolchain)

    # eclipse specific flag
    env.Append(CXXFLAGS=' -fmessage-length=0')

    # c++0x flag
    env.Append(CXXFLAGS=' -std=c++0x')

    # support __float128
#    env.Append(CXXFLAGS=' -D__STRICT_ANSI__')

    # be pedantic to avoid stupid programming errors
    env.Append(CXXFLAGS=' -pedantic')

    # speedup compilation - remove this when compiler slows down or segfaults by running out of memory
    env.Append(CXXFLAGS=' -pipe')

    env.Replace(CXX = 'clang++')
    env.Replace(LINK = 'clang++')

    if p.fortran_source == 'enable':
        env.Replace(F90='gfortran')
        env.Append(F90FLAGS=' -cpp')
        env.Append(LIBS=['gfortran'])



if p.mode in ['debug', 'debug_thread', 'debug_leak']:
    env.Append(CXXFLAGS=' -DSWEET_DEBUG=1')

    if compiler_to_use == 'gnu':
        env.Append(CXXFLAGS='-O0 -g3 -Wall')

        # integer overflow check
        env.Append(CXXFLAGS=' -ftrapv')

    elif compiler_to_use == 'llvm':
        env.Append(CXXFLAGS='-O0 -g3 -Wall')

    elif compiler_to_use == 'intel':
        env.Append(CXXFLAGS='-O0 -g -traceback')
#        env.Append(CXXFLAGS=' -fp-trap=common')

    elif compiler_to_use == 'pgi':
        env.Append(CXXFLAGS='-O0 -g -traceback')


    if p.fortran_source == 'enable':
        if compiler_to_use == 'gnu':
            env.Append(F90FLAGS='-g -O0')
        elif compiler_to_use == 'intel':
            env.Append(F90FLAGS='-g -O0 -traceback')


elif p.mode == 'release':
    env.Append(CXXFLAGS=' -DSWEET_DEBUG=0')

    # deactivate assertion calls
    env.Append(CXXFLAGS=' -DNDEBUG=1')

    if compiler_to_use == 'gnu':
        env.Append(CXXFLAGS=' -O3 -mtune=native')

        # Ensure vectorization
        env.Append(CXXFLAGS=' -ftree-vectorize')

        # Let the compiler know about no aliasing of function arguments
        env.Append(CXXFLAGS=' -fstrict-aliasing')

    elif compiler_to_use == 'llvm':
        env.Append(CXXFLAGS=' -O3 -mtune=native')

    elif compiler_to_use == 'intel':
        env.Append(CXXFLAGS=' -O2 -fno-alias')

        if p.mic != 'enable':
            env.Append(CXXFLAGS=' -xHost')

    elif compiler_to_use == 'pgi':
        env.Append(CXXFLAGS='-O3 -fast -Mipa=fast,inline -Msmartalloc')

    if p.fortran_source == 'enable':
        if compiler_to_use == 'gnu':
            env.Append(F90FLAGS=' -O2')
        elif compiler_to_use == 'intel':
            env.Append(F90FLAGS=' -O2')

if p.sanitize != '':
    env.Append(CXXFLAGS=' -fsanitize='+p.sanitize)


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
    env.ParseConfig("pkg-config --libs --cflags freetype2")

else:
    env.Append(CXXFLAGS=' -DSWEET_GUI=0')

  

if p.fortran_source == 'enable':
    env.Append(CXXFLAGS=' -DSWEET_FORTRAN=1')
else:
    env.Append(CXXFLAGS=' -DSWEET_FORTRAN=0')


if exec_command('uname -s') == "Darwin":
    env.Append(CXXFLAGS=' -DMEMBLOCKALLOC_ENABLE_NUMA_ALLOC=0')
else:
    env.Append(LIBS=['numa'])


if p.threading == 'omp':
    env.Append(CXXFLAGS=['-fopenmp'])
    env.Append(LINKFLAGS=['-fopenmp'])

    env.Append(CXXFLAGS=' -DSWEET_THREADING_SPACE=1')
else:
    env.Append(CXXFLAGS=' -DSWEET_THREADING_SPACE=0')



if p.eigen == 'enable':
    env.Append(CXXFLAGS=['-Ilocal_software/local/include/eigen3'])
    env.Append(CXXFLAGS=['-DSWEET_EIGEN=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_EIGEN=0'])






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

    if compiler_to_use == 'gnu':
        env.Append(LIBS=['gfortran'])

    p.libfft = 'enable'



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
        env.Append(LIBS=['fftw3'])

        if p.threading == 'omp' or p.rexi_thread_parallel_sum == 'enable':
            env.Append(LIBS=['fftw3_omp'])

else:
    env.Append(CXXFLAGS = ' -DSWEET_USE_LIBFFT=0')

if p.mic == 'enable':
    env.Append(CXXFLAGS=['-mmic'])
    env.Append(LINKFLAGS=['-mmic'])


#if p.rexi_thread_parallel_sum == 'enable' and p.threading == 'omp':
#    raise Exception('ERROR: "REXI Parallel Sum" and "Threading" is both activated')



#
# If SWEET_THREADING_TIME_REXI is activated, the REXI sum is computed
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
    #if SWEET_THREADING_SPACE
    #
    env.Append(CXXFLAGS=['-fopenmp'])

    #
    # Compile flag is set in sconscript
    #

    # Activate precompiler flag
    env.Append(CXXFLAGS=' -DSWEET_THREADING_TIME_REXI=1')
else:
    env.Append(CXXFLAGS=' -DSWEET_THREADING_TIME_REXI=0')


if p.benchmark_timings == 'enable':
    env.Append(CXXFLAGS=' -DSWEET_BENCHMARK_TIMINGS=1')
else:
    env.Append(CXXFLAGS=' -DSWEET_BENCHMARK_TIMINGS=0')

if p.rexi_timings_additional_barriers == 'enable':
    env.Append(CXXFLAGS=' -DSWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS=1')
else:
    env.Append(CXXFLAGS=' -DSWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS=0')

if p.rexi_allreduce == 'enable':
    env.Append(CXXFLAGS=' -DSWEET_REXI_ALLREDUCE=1')
else:
    env.Append(CXXFLAGS=' -DSWEET_REXI_ALLREDUCE=0')


if p.threading == 'omp' or p.rexi_thread_parallel_sum == 'enable':
    env.Append(CXXFLAGS=' -DSWEET_THREADING=1')
else:
    env.Append(CXXFLAGS=' -DSWEET_THREADING=0')


if p.debug_symbols == 'enable':
    env.Append(CXXFLAGS=['-g'])
    env.Append(LINKFLAGS=['-g'])

    if compiler_to_use == 'intel':
        env.Append(CXXFLAGS=['-shared-intel -shared-libgcc -debug inline-debug-info'])
        env.Append(LINKFLAGS=['-shared-intel  -shared-libgcc -debug inline-debug-info'])

if compiler_to_use != p.compiler:
    print("*"*80)
    print("WARNING: Using compiler "+str(compiler_to_use))
    print("*"*80)

exec_name = p.getProgramName()


#
# BUILD directory
#
build_dir='/tmp/'
user = os.environ.get('USERNAME')
if user != None:
    build_dir += user.replace(' ', '')
    build_dir += '/'

build_dir += 'scons_build_'+exec_name+'/'

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

    env.SConscript('sconscript.py', variant_dir=build_dir, duplicate=0, exports=['env', 'p'])

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
