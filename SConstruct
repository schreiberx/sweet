import os
import subprocess
import multiprocessing

from mule.JobCompileOptions import *

#
# Setup parallel compilation
#
num_cpu = int(os.environ.get('MULE_COMPILE_NUM_JOB_LIMITATION', -1))
if num_cpu == -1:
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

jg = JobCompileOptions()



###################################################################
# Determine compiler to use
###################################################################


compiler_type_cxx = None
compiler_type_f90 = None

if 'MULE_CXX_COMPILER' in os.environ:
    compiler_type_cxx = os.environ['MULE_CXX_COMPILER']

else:
    print("MULE_CXX_COMPILER env not found, trying to autodetect compiler")

    if jg.sweet_mpi == 'enable':
        raise Exception("Please specify MULE_CXX_COMPILER with MPI to ensure no compile problems")

    if 'CXX' in os.environ:
        cxx = os.environ['CXX']
        if cxx.startswith("g++"):
            compiler_type_cxx = 'gcc'
        elif cxx.startswith("clang"):
            compiler_type_cxx = 'llvm'


if 'MULE_F90_COMPILER' in os.environ:
    compiler_type_f90 = os.environ['MULE_F90_COMPILER']

else:
    print("MULE_F90_COMPILER env not found, trying to autodetect compiler")

    if jg.sweet_mpi == 'enable':
        raise Exception("Please specify MULE_F90_COMPILER with MPI to ensure no compile problems")

    if 'F90' in os.environ:
        f90 = os.environ['F90']
        if f90.startswith("gfortran"):
            compiler_type_f90 = 'gcc'
        elif f90.startswith("flang"):
            compiler_type_f90 = 'llvm'

if compiler_type_cxx is None:
    raise Exception("Failed to automatically detect C++ compiler")

print(f"Using CXX compiler '{compiler_type_cxx}'")
print(f"Using F90 compiler '{compiler_type_f90}'")



###################################################################
# fix LD LIB PATH
###################################################################
if 'LD_LIBRARY_PATH' in os.environ:
    env.Append(LIBPATH=os.environ['LD_LIBRARY_PATH'].split(':'))



###########################################
# Compile options
###########################################


jg.sconsProcessCommandlineOptions()

# Cleanup options
jg.postprocessOptions()



if jg.libxml == 'enable':
    env.ParseConfig("xml2-config --cflags --libs")


env.Append(CXXFLAGS=['-DSWEET_SIMD_ENABLE='+('1' if jg.simd=='enable' else '0')])
env.Append(CXXFLAGS=['-DCONFIG_ENABLE_LIBXML='+('1' if jg.libxml=='enable' else '0')])



if jg.plane_spectral_space == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_USE_PLANE_SPECTRAL_SPACE=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_USE_PLANE_SPECTRAL_SPACE=0'])



if jg.plane_spectral_dealiasing == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_USE_PLANE_SPECTRAL_DEALIASING=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_USE_PLANE_SPECTRAL_DEALIASING=0'])



if jg.sphere_spectral_space == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_USE_SPHERE_SPECTRAL_SPACE=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_USE_SPHERE_SPECTRAL_SPACE=0'])


if jg.sphere_spectral_dealiasing == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_USE_SPHERE_SPECTRAL_DEALIASING=1'])
else:
    if jg.sphere_spectral_space == 'enable':
        raise Exception("No anti-aliasing on sphere as compile option supported, please simply use command line options to specify lower physical resolution!")

    env.Append(CXXFLAGS=['-DSWEET_USE_SPHERE_SPECTRAL_DEALIASING=0'])

#
# Add other flags, potentially from shacks
#
jg.sconsAddFlags(env)


#
# Override compiler settings from environment variable if one of these variables is set
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

    if jg.sweet_mpi == 'enable':
        mi = 'MULE_MPI'+i
        if mi in env['ENV']:
            if 'FLAGS' in i or 'LIBS' in i:
                print("INFO: Appending to "+i+"+= "+env['ENV'][mi])
                env.Append(**{i: [env['ENV'][mi]]})
            else:
                print("INFO: Using MULE_MPI* environment variable to set "+i+"="+env['ENV'][mi])
                env[i] = env['ENV'][mi]


    else:
        if i in env['ENV']:
            if 'FLAGS' in i or 'LIBS' in i:
                print("INFO: Appending to "+i+"+= "+env['ENV'][i])
                env.Append(**{i: [env['ENV'][i]]})

            else:
                print("INFO: Overriding environment variable "+i+"="+env['ENV'][i])
                env[i] = env['ENV'][i]



#
# MPI
#
if jg.sweet_mpi == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_MPI=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_MPI=0'])


if jg.libpfasst == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_LIBPFASST=1'])

    env.Append(LIBS=['libpfasst'])

    # enable MPI per default for libpfasst
    jg.sweet_mpi = 'enable'

else:
    env.Append(CXXFLAGS=['-DSWEET_LIBPFASST=0'])



env.Append(LIBS=['m'])


if compiler_type_cxx == 'gcc':

    # eclipse specific flag
    env.Append(CXXFLAGS=['-fmessage-length=0'])

    # c++0x flag
    env.Append(CXXFLAGS=['-std=c++0x'])


    if jg.sweet_mpi == 'enable':
        # GNU compiler needs special treatment!
        # Linking with Fortran MPI requires
        # for OpenMPI: -lmpi_mpifh
        # for MPICH: -lmpifort

        output = exec_command(env['CC']+' -v')
        if 'MPICH' in output:
            if jg.fortran_source == 'enable':
                env.Append(LIBS=['mpifort'])
        else:
            # Assume OpenMPI
            if jg.fortran_source == 'enable':
                env.Append(LIBS=['mpi_cxx'])


elif compiler_type_cxx == 'intel':
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

    env.Append(LINKFLAGS=['-shared-intel'])
    env.Append(LINKFLAGS=['-shared-libgcc'])
    env.Append(LINKFLAGS=['-debug inline-debug-info'])

    # eclipse specific flag
    env.Append(CXXFLAGS=['-fmessage-length=0'])

    # c++0x flag
    env.Append(CXXFLAGS=['-std=c++0x'])

    # output more warnings
    env.Append(CXXFLAGS=['-w1'])

    if jg.sweet_mpi == 'enable':
        if jg.threading != 'off':
            env.Append(CXXFLAGS=['-mt_mpi'])
            env.Append(LINKFLAGS=['-mt_mpi'])


elif compiler_type_cxx == 'llvm':
    reqversion = [9,0]
    if jg.threading == 'omp':
        reqversion = [9,0]
    version_line = exec_command(f"{os.environ['CXX']} --version").splitlines()[0]

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

    # eclipse specific flag
    env.Append(CXXFLAGS=['-fmessage-length=0'])

    # c++0x flag
    env.Append(CXXFLAGS=['-std=c++0x'])

    # support __float128
#    env.Append(CXXFLAGS='-D__STRICT_ANSI__')

    # be pedantic to avoid stupid programming errors
    env.Append(CXXFLAGS=['-pedantic'])

    # speedup compilation - remove this when compiler slows down or segfaults by running out of memory
    env.Append(CXXFLAGS=['-pipe'])

else:
    raise Exception("Unsupported compiler")



"""
Different modes
"""

if jg.mode in ['debug', 'debug_thread', 'debug_leak']:
    env.Append(CXXFLAGS=['-DSWEET_DEBUG=1'])

    if compiler_type_cxx == 'gcc':
        env.Append(CXXFLAGS=["-O0", "-g3", "-Wall"])

        # integer overflow check
        env.Append(CXXFLAGS=['-ftrapv'])

        #env.Append(CXXFLAGS=['-fsanitize=address', '-fsanitize=undefined', '-fno-sanitize-recover=all', '-fsanitize=float-divide-by-zero', '-fsanitize=float-cast-overflow', '-fno-sanitize=null', '-fno-sanitize=alignment'])

    elif compiler_type_cxx == 'llvm':
        env.Append(CXXFLAGS=["-O0", "-g3", "-Wall"])

        # Memory sanitizer
        #env.Append(CXXFLAGS=['-fsanitize=memory'])
        #env.Append(CXXFLAGS=['-fsanitize=address', '-fsanitize=undefined', '-fno-sanitize-recover=all', '-fsanitize=float-divide-by-zero', '-fsanitize=float-cast-overflow', '-fno-sanitize=null', '-fno-sanitize=alignment'])

    elif compiler_type_cxx == 'intel':
        env.Append(CXXFLAGS=["-O0", "-g", "-traceback"])
#        env.Append(CXXFLAGS='-fp-trap=common')


    if jg.fortran_source == 'enable':
        if compiler_type_cxx == 'gcc':
            env.Append(F90FLAGS=["-g", "-O0"])
        elif compiler_type_cxx == 'intel':
            env.Append(F90FLAGS=["-g", "-O0", "-traceback"])


elif jg.mode == 'release':
    env.Append(CXXFLAGS=['-DSWEET_DEBUG=0'])

    # deactivate assertion calls
    env.Append(CXXFLAGS=['-DNDEBUG=1'])

    if compiler_type_cxx == 'gcc':
        env.Append(CXXFLAGS=["-O3", "-mtune=native", "-march=native"])

        # Ensure vectorization
        env.Append(CXXFLAGS=['-ftree-vectorize'])

        # Let the compiler know about no aliasing of function arguments
        env.Append(CXXFLAGS=['-fstrict-aliasing'])

    elif compiler_type_cxx == 'llvm':
        env.Append(CXXFLAGS=["-O3", "-mtune=native", "-march=native"])

    elif compiler_type_cxx == 'intel':
        env.Append(CXXFLAGS=["-O2", "-fno-alias"])

        if jg.mic != 'enable':
            env.Append(CXXFLAGS=['-xHost'])

    if jg.fortran_source == 'enable':
        if compiler_type_cxx == 'gcc':
            env.Append(F90FLAGS=["-O2"])
        elif compiler_type_cxx == 'intel':
            env.Append(F90FLAGS=["-O2"])


if jg.quadmath == 'enable':
    env.Append(CXXFLAGS=["-DSWEET_QUADMATH=1"])
    env.Append(LIBS=['quadmath'])
else:
    env.Append(CXXFLAGS=["-DSWEET_QUADMATH=0"])



if jg.gui == 'enable':
    # compile flags
    env.Append(CXXFLAGS=['-I'+os.environ['HOME']+'/local/include'])
    env.Append(CXXFLAGS=['-DSWEET_GUI=1'])


    if exec_command('uname -s') == "Darwin":
        # ASSUME MACOSX SYSTEM
        env.Append(LINKFLAGS=['-framework OpenGL'])
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
    env.Append(CXXFLAGS=['-DSWEET_GUI=0'])



# Add LAPACK libraries
if jg.lapack == 'enable':
    env.Append(LIBS=['lapack'])
    env.Append(LIBS=['blas'])
    env.Append(CXXFLAGS=['-DSWEET_LAPACK=1'])

    assert jg.fortran_source == 'enable', "LAPACK enabled, Fortran source needs to be enabled for lapack (--fortran-source=enbale) "

else:
    env.Append(CXXFLAGS=['-DSWEET_LAPACK=0'])


  

if jg.fortran_source == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_FORTRAN=1'])

    if compiler_type_cxx == 'gcc':
        env.Append(LIBS=['gfortran'])

else:
    env.Append(CXXFLAGS=['-DSWEET_FORTRAN=0'])


if exec_command("uname -s").replace("\n", "") == "Darwin":
    env.Append(CXXFLAGS=['-DMEMBLOCKALLOC_ENABLE_NUMA_ALLOC=0'])
else:
    env.Append(LIBS=['numa'])


if jg.threading == 'omp':
    env.Append(CXXFLAGS=['-fopenmp'])
    env.Append(LINKFLAGS=['-fopenmp'])

    env.Append(CXXFLAGS=['-DSWEET_THREADING_SPACE=1'])
elif jg.threading == 'off':
    env.Append(CXXFLAGS=['-DSWEET_THREADING_SPACE=0'])
else:
    raise Exception("Internal error")



if jg.eigen == 'enable':
    env.Append(CXXFLAGS=['-Ilocal_software/local/include/eigen3'])
    env.Append(CXXFLAGS=['-DSWEET_EIGEN=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_EIGEN=0'])



if jg.libsph == 'enable':
    if jg.threading == 'omp':
        env.Append(LIBS=['shtns_omp'])
    else:
        env.Append(LIBS=['shtns'])

    jg.libfft = 'enable'



if jg.libfft == 'enable':

    env.Append(CXXFLAGS='-DSWEET_USE_LIBFFT=1')

    if jg.mkl == 'enable':
        print("INFO: Using Intel MKL instead of FFTW");

        if jg.threading != 'omp':
            env.Append(CXXFLAGS=['-mkl=sequential'])
            env.Append(LINKFLAGS=['-mkl=sequential'])
        else:
            env.Append(CXXFLAGS=['-mkl=parallel'])
            env.Append(LINKFLAGS=['-mkl=parallel'])

    else:
        env.Append(LIBS=['fftw3'])

        if jg.threading == 'omp' or jg.rexi_thread_parallel_sum == 'enable':
            env.Append(LIBS=['fftw3_omp'])

else:
    env.Append(CXXFLAGS=['-DSWEET_USE_LIBFFT=0'])


#
# If SWEET_THREADING_TIME_REXI is activated, the REXI sum is computed
# with parallel for over the sum terms
#
if jg.rexi_thread_parallel_sum == 'enable':
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
    env.Append(CXXFLAGS=['-DSWEET_THREADING_TIME_REXI=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_THREADING_TIME_REXI=0'])

if jg.parallel_sdc_par_model == 'omp':
    # Same for gcc/icpc
    env.Append(LINKFLAGS=['-fopenmp'])
    env.Append(CXXFLAGS=['-fopenmp'])

    # Activate precompiler flag
    env.Append(CXXFLAGS=['-DSWEET_PARALLEL_SDC_OMP_MODEL=1'])


if jg.benchmark_timings == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_BENCHMARK_TIMINGS=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_BENCHMARK_TIMINGS=0'])

if jg.rexi_timings_additional_barriers == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_REXI_TIMINGS_ADDITIONAL_BARRIERS=0'])

if jg.rexi_allreduce == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_REXI_ALLREDUCE=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_REXI_ALLREDUCE=0'])


if jg.threading == 'omp' or jg.rexi_thread_parallel_sum == 'enable':
    env.Append(CXXFLAGS=['-DSWEET_THREADING=1'])
else:
    env.Append(CXXFLAGS=['-DSWEET_THREADING=0'])


if jg.debug_symbols == 'enable':
    env.Append(CXXFLAGS=['-g'])
    env.Append(LINKFLAGS=['-g'])

    if compiler_type_cxx == 'intel':
        env.Append(CXXFLAGS=['-shared-intel', '-shared-libgcc', '-debug',  'inline-debug-info'])
        env.Append(LINKFLAGS=['-shared-intel', '-shared-libgcc', '-debug',  'inline-debug-info'])


"""
Fortran related support at the very end here
"""
if compiler_type_cxx == 'gcc':
    if jg.fortran_source == 'enable':
        env.Append(F90FLAGS=['-cpp'])
        env.Append(LIBS=['gfortran'])

elif compiler_type_cxx == 'intel':
    if jg.fortran_source == 'enable':
        env.Append(LIBS=['gfortran'])
        env.Append(LIBS=['ifcore'])
        env.Append(F90FLAGS=['-fpp'])

elif compiler_type_cxx == 'llvm':
    if jg.fortran_source == 'enable':
        env.Append(F90FLAGS=['-cpp'])
        env.Append(LIBS=['gfortran'])



exec_name = jg.getProgramExec()


#
# BUILD directory
#
build_dir='/tmp/'
user = os.environ.get('USERNAME')
if user != None:
    build_dir += user.replace(' ', '')
    build_dir += '/'

build_dir += 'scons_build_'+exec_name+'/'

if jg.libpfasst == 'enable':
    env.Append(F90FLAGS = ['-Ilocal_software/local/include/libpfasst'])

#
# USE build directory for Fortran module output
#
env.Append(F90FLAGS = ['-J'+build_dir])


# Leave this deactivated!
# env.Append(CPPPATH = ['/usr/local/include'])




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


if jg.program_name != 'DUMMY':

    env.SConscript('sconscript.py', variant_dir=build_dir, duplicate=0, exports=['env', 'jg'])

    print('')
    print('            Program: '+jg.program_name)
    print('Building executable: '+exec_name)
    print('')

    obj_files = []
    for i in env.src_files:
        for a in i:
            if str(a).split(".")[-1] == "o":
                obj_files.append(str(a))

    env.Program('build/'+exec_name, obj_files)

