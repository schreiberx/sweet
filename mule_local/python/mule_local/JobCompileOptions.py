#! /usr/bin/env python3
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
from mule.InfoError import *

import subprocess


__all__ = ['JobCompileOptions']


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


class JobCompileOptions(InfoError):


    def __init__(self, dummy_init = False):
        self.init_phase = True

        InfoError.__init__(self, "JobCompileOptions")

        # Program or unit test
        self.program = ""
        self.program_name = ""
        self.unit_test = ""

        # Compile options
        self.mode = 'release'

        self.debug_symbols = 'disable'
        self.simd = 'enable'

        self.fortran_source = 'disable'

        self.lapack = 'enable'

        self.program_binary_name = ''

        # Parallelization
        self.sweet_mpi = 'disable'
        self.threading = 'omp'
        self.rexi_thread_parallel_sum = 'disable'

        self.benchmark_timings = 'disable'

        # Additional barriers to overcome issues of turbo boost
        self.rexi_timings_additional_barriers = 'disable'

        # Use reduce all instead of reduce to root rank
        self.rexi_allreduce = 'disable'


        # Program / Unit test
        self.program = ''
        self.unit_test = ''

        # Parareal
        self.parareal = 'none'
        self.parareal_scalar = 'disable'
        self.parareal_plane = 'disable'
        self.parareal_sphere = 'disable'
        self.parareal_plane_swe = 'disable'
        self.parareal_plane_burgers = 'disable'

        # XBraid
        self.xbraid = 'none'
        self.xbraid_scalar = 'disable'
        self.xbraid_plane = 'disable'
        self.xbraid_sphere = 'disable'
        self.xbraid_plane_swe = 'disable'
        self.xbraid_plane_burgers = 'disable'

        #LibPFASST
        self.libpfasst = 'disable'

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

        self.init_phase = False



    def __setattr__(self, name, value):
        if name != 'init_phase':
            if not self.init_phase:
                if not name in self.__dict__:
                    raise Exception("Attribute '"+name+"' does not exist!")

        self.__dict__[name] = value



    def getSConsParams(self):
        retval = ''
        retval += ' --mode='+self.mode
        retval += ' --debug-symbols='+("enable" if self.debug_symbols else "disable")
        retval += ' --simd='+self.simd

        retval += ' --fortran-source='+self.fortran_source

        retval += ' --lapack='+self.lapack

        retval += ' --program-binary-name='+self.program_binary_name

        # Parallelization
        retval += ' --sweet-mpi='+self.sweet_mpi
        retval += ' --threading='+self.threading
        retval += ' --rexi-thread-parallel-sum='+self.rexi_thread_parallel_sum
        retval += ' --benchmark-timings='+self.benchmark_timings
        retval += ' --rexi-timings-additional-barriers='+self.rexi_timings_additional_barriers
        retval += ' --rexi-allreduce='+self.rexi_allreduce

        # Program / Unit test
        if self.program != '':
            retval += ' --program='+self.program
        if self.unit_test != '':
            retval += ' --unit-test='+self.unit_test

        # Parareal
        retval += ' --parareal='+self.parareal
        retval += ' --parareal-scalar='+self.parareal_scalar
        retval += ' --parareal-plane='+self.parareal_plane
        retval += ' --parareal-sphere='+self.parareal_sphere
        retval += ' --parareal-plane-swe='+self.parareal_plane_swe
        retval += ' --parareal-plane-burgers='+self.parareal_plane_burgers

        # XBraid
        retval += ' --xbraid='+self.xbraid
        retval += ' --xbraid-scalar='+self.xbraid_scalar
        retval += ' --xbraid-plane='+self.xbraid_plane
        retval += ' --xbraid-sphere='+self.xbraid_sphere
        retval += ' --xbraid-plane-swe='+self.xbraid_plane_swe
        retval += ' --xbraid-plane-burgers='+self.xbraid_plane_burgers


        # LibPFASST
        retval += ' --libpfasst='+self.libpfasst

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

#        scons.AddOption(    '--xml-config',
#                dest='xml_config',
#                type='string',
#                default='',
#                help='xml configuration file with compile options'
#        )
#        env['xml_config'] = scons.GetOption('xml_config')
#
#        if env['xml_config'] != '':
#            self.env = CompileXMLOptions.load(env['xml_config'], env)



        scons.AddOption(      '--mode',
                dest='mode',
                type='choice',
                choices=['debug', 'release'],
                default='release',
                help='specify compilation mode to use: debug, release [default: %default]'
        )
        self.mode = scons.GetOption('mode')

        scons.AddOption(    '--simd',
                dest='simd',
                type='choice',
                choices=['enable', 'disable'],
                default='enable',
                help="Use SIMD for operations such as folding [default: %default]"
        )
        self.simd = scons.GetOption('simd')


        scons.AddOption(    '--eigen',
                dest='eigen',
                type='choice',
                choices=['enable', 'disable'],
                default='disable',
                help="Activate utilization of Eigen library [default: %default]"
        )
        self.eigen = scons.GetOption('eigen')


        scons.AddOption(    '--libpfasst',
                dest='libpfasst',
                type='choice',
                choices=['enable', 'disable'],
                default='disable',
                help="Activate utilization of libPFASST (FOortran version) [default: %default]"
        )
        self.libpfasst = scons.GetOption('libpfasst')



        scons.AddOption(    '--debug-symbols',
                dest='debug_symbols',
                type='choice',
                choices=['enable', 'disable'],
                default='enable',
                help="Create binary with debug symbols [default: %default]"
        )
        self.debug_symbols = scons.GetOption('debug_symbols')



        scons.AddOption(    '--plane-spectral-space',
                dest='plane_spectral_space',
                type='choice',
                choices=['enable', 'disable'],
                default='enable',
                help="Activate spectral space for data on the plane (2D FFT) [default: %default]"
        )
        self.plane_spectral_space = scons.GetOption('plane_spectral_space')

        scons.AddOption(    '--sphere-spectral-space',
                dest='sphere_spectral_space',
                type='choice',
                choices=['enable', 'disable'],
                default='disable',
                help="Activate spectral space for data on the sphere (Spherical Harmonics) [default: %default]"
        )
        self.sphere_spectral_space = scons.GetOption('sphere_spectral_space')


        scons.AddOption(    '--fortran-source',
                dest='fortran_source',
                type='choice',
                choices=['enable', 'disable'],
                default='disable',
                help="Activate linking with Fortran source [default: %default]"
        )
        self.fortran_source = scons.GetOption('fortran_source')


        scons.AddOption(    '--lapack',
                dest='lapack',
                type='choice',
                choices=['enable', 'disable'],
                default='enable',
                help="Enable lapack [default: %default]"
        )
        self.lapack = scons.GetOption('lapack')



        scons.AddOption(    '--libfft',
                dest='libfft',
                type='choice',
                choices=['enable', 'disable'],
                default='disable',
                help="Enable compiling and linking with FFT library [default: %default]"
        )
        self.libfft = scons.GetOption('libfft')


        scons.AddOption(    '--libsph',
                dest='libsph',
                type='choice',
                choices=['enable', 'disable'],
                default='disable',
                help="Enable compiling and linking with SPH library [default: %default]"
        )
        self.libsph = scons.GetOption('libsph')


        scons.AddOption(    '--mkl',
                dest='mkl',
                type='choice',
                choices=['enable', 'disable'],
                default='disable',
                help="Enable Intel MKL [default: %default]"
        )
        self.mkl = scons.GetOption('mkl')


        #
        # LIB XML
        #
        scons.AddOption(    '--libxml',
                dest='libxml',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='Compile with libXML related code: enable, disable [default: %default]'
        )
        self.libxml = scons.GetOption('libxml')


        scons.AddOption(    '--plane-spectral-dealiasing',
                dest='plane_spectral_dealiasing',
                type='choice',
                choices=['enable','disable'],
                default='enable',
                help='spectral dealiasing (3N/2-1 rule): enable, disable [default: %default]'
        )
        self.plane_spectral_dealiasing = scons.GetOption('plane_spectral_dealiasing')

        scons.AddOption(    '--sphere-spectral-dealiasing',
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

        scons.AddOption(    '--quadmath',
                dest='quadmath',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='quadmath: enable, disable [default: %default]'
        )
        self.quadmath = scons.GetOption('quadmath')



        scons.AddOption(    '--gui',
                dest='gui',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='gui: enable, disable [default: %default]'
        )
        self.gui = scons.GetOption('gui')




        scons.AddOption(    '--rexi-thread-parallel-sum',
                dest='rexi_thread_parallel_sum',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='Use a par for loop over the sum in REXI: (enable, disable) [default: %default]\n\tWARNING: This also disables the parallelization-in-space with OpenMP'
        )
        self.rexi_thread_parallel_sum = scons.GetOption('rexi_thread_parallel_sum')


        scons.AddOption(    '--benchmark-timings',
                dest='benchmark_timings',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='REXI timings: enable, disable [default: %default]'
        )
        self.benchmark_timings = scons.GetOption('benchmark_timings')

        scons.AddOption(    '--rexi-timings-additional-barriers',
                dest='rexi_timings_additional_barriers',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='REXI timings with additional barriers: enable, disable [default: %default]\nThis is helpful for improved measurements if TurboBoost is activated'
        )
        self.rexi_timings_additional_barriers = scons.GetOption('rexi_timings_additional_barriers')

        scons.AddOption(    '--rexi-allreduce',
                dest='rexi_allreduce',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='For REXI, use MPI_Allreduce operations instead of MPI_Reduce: enable, disable [default: %default]'
        )
        self.rexi_allreduce = scons.GetOption('rexi_allreduce')


        scons.AddOption(    '--sweet-mpi',
                dest='sweet_mpi',
                type='choice',
                choices=['enable','disable'],
                default='disable',
                help='Enable MPI commands e.g. for REXI sum: enable, disable [default: %default]'
        )
        self.sweet_mpi = scons.GetOption('sweet_mpi')



        scons.AddOption(    '--parareal',
                dest='parareal',
                type='choice',
                choices=['none', 'serial','mpi'],
                default='none',
                help='Enable Parareal (none, serial, mpi) [default: %default]\nOnly works, if Parareal is supported by the simulation'
        )
        self.parareal = scons.GetOption('parareal')

        scons.AddOption(    '--parareal-scalar',
                dest='parareal_scalar',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal for scalar problems (enable, disable) [default: %default]'
        )
        self.parareal_scalar = scons.GetOption('parareal_scalar')

        scons.AddOption(    '--parareal-plane',
                dest='parareal_plane',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal on the plane (enable, disable) [default: %default]'
        )
        self.parareal_plane = scons.GetOption('parareal_plane')

        scons.AddOption(    '--parareal-sphere',
                dest='parareal_sphere',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal on the sphere (enable, disable) [default: %default]'
        )
        self.parareal_sphere = scons.GetOption('parareal_sphere')

        scons.AddOption(    '--parareal-plane-swe',
                dest='parareal_plane_swe',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal for SWE on the plane (enable, disable) [default: %default]'
        )
        self.parareal_plane_swe = scons.GetOption('parareal_plane_swe')

        scons.AddOption(    '--parareal-plane-burgers',
                dest='parareal_plane_burgers',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal for Burgers on the plane (enable, disable) [default: %default]'
        )
        self.parareal_plane_burgers = scons.GetOption('parareal_plane_burgers')

        scons.AddOption(    '--xbraid',
                dest='xbraid',
                type='choice',
                choices=['none', 'mpi'],
                default='none',
                help='Enable XBBraid (none,  mpi) [default: %default]\nOnly works, if XBraid is supported by the simulation'
        )
        self.xbraid = scons.GetOption('xbraid')

        scons.AddOption(    '--xbraid-scalar',
                dest='xbraid_scalar',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for scalar problems (enable, disable) [default: %default]'
        )
        self.xbraid_scalar = scons.GetOption('xbraid_scalar')

        scons.AddOption(    '--xbraid-plane',
                dest='xbraid_plane',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid on the plane (enable, disable) [default: %default]'
        )
        self.xbraid_plane = scons.GetOption('xbraid_plane')

        scons.AddOption(    '--xbraid-sphere',
                dest='xbraid_sphere',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid on the sphere (enable, disable) [default: %default]'
        )
        self.xbraid_sphere = scons.GetOption('xbraid_sphere')

        scons.AddOption(    '--xbraid-plane-swe',
                dest='xbraid_plane_swe',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for SWE on the plane (enable, disable) [default: %default]'
        )
        self.xbraid_plane_swe = scons.GetOption('xbraid_plane_swe')

        scons.AddOption(    '--xbraid-plane-burgers',
                dest='xbraid_plane_burgers',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for Burgers on the plane (enable, disable) [default: %default]'
        )
        self.xbraid_plane_burgers = scons.GetOption('xbraid_plane_burgers')



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


        if not self.parareal == 'none':
            if self.program == "parareal_ode":
                self.parareal_scalar = 'enable';
            elif self.program == 'swe_plane' or self.program == 'burgers':
                self.parareal_plane = 'enable';
                if self.program == 'swe_plane':
                    self.parareal_plane_swe = 'enable';
                elif self.program == 'burgers':
                    self.parareal_plane_burgers = 'enable';
            elif self.program == 'swe_sphere':
                self.parareal_sphere = 'enable';

        if not self.xbraid == 'none':
            if self.program == "parareal_ode":
                self.xbraid_scalar = 'enable';
            elif self.program == 'swe_plane' or self.program == 'burgers':
                self.xbraid_plane = 'enable';
                if self.program == 'swe_plane':
                    self.xbraid_plane_swe = 'enable';
                elif self.program == 'burgers':
                    self.xbraid_plane_burgers = 'enable';
            elif self.program == 'swe_sphere':
                self.xbraid_sphere = 'enable';


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
        scons.AddOption(    '--threading',
                dest='threading',
                type='choice',
                choices=threading_constraints,
                default='omp',
                help='Threading to use '+' / '.join(threading_constraints)+', default: off'
        )
        self.threading = scons.GetOption('threading')


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

        retval = self.program_name
        s = self.getUniqueID([])

        if s != '':
            retval += '_'+s

        return retval


    #
    # Return a dictionary with further compile options specified in the
    # headers of the main program file
    #
    def get_program_specific_options(self):

        mainsrcadddir = ''

        if self.program != '':
            mainsrcadddir = 'src/programs/'+self.program
        elif self.unit_test != '':
            mainsrcadddir = 'src/unit_tests/'+self.unit_test
        else:
            return None


        #
        # Add main source file
        #
        main_src = mainsrcadddir+'.cpp'

        #
        # Check for additional source files and directories which should be added
        # This can be specified in the main program/test source file via e.g.
        #
        # MULE_COMPILE_DIRS: [file1] [dir2] [file2] [file3]
        #
        # The files and directories need to be specified to the software's root directory
        #

        sw_root = os.environ['MULE_SOFTWARE_ROOT']+'/'
        f = open(sw_root+'/'+main_src, 'r')
        lines = f.readlines()

        fad_dict = {
            'compile_files_and_dirs': [],
            'scons_options': [],
        }

        tags_ = [
                    ['compile_files_and_dirs', 'MULE_COMPILE_FILES_AND_DIRS: '],
                    ['scons_options', 'MULE_SCONS_OPTIONS: '],
                ]

        # Add main program file to compile files
        fad_dict['compile_files_and_dirs'] += [main_src]

        for l in lines:
            for tags in tags_:
                tag_id = tags[0]
                tag = tags[1]
                if tag in l:
                    fad = l[l.find(tag)+len(tag):]
                    fad = fad.replace('\r', '')
                    fad = fad.replace('\n', '')

                    fad_dict[tag_id] += fad.split(' ')

        return fad_dict
        


    def process_scons_options(self, options_):
        for option in options_:
            name, value = option.split('=', 1)

            if name[0:2] != '--':
                raise Exception("Option doesn't start with '--': "+name)

            name = name[2:]

            if name == 'sphere-spectral-space':
                self.sphere_spectral_space = value
                continue

            if name == 'plane-spectral-space':
                self.plane_spectral_space = value
                continue

            if name == 'fortran-source':
                self.fortran_source = value
                continue

            raise Exception("TODO: Process option '"+name+"'. Maybe you need to add it here!")


    def getProgramPath(self, ignore_errors = False):
        return os.environ['MULE_SOFTWARE_ROOT']+'/build/'+self.getProgramName(ignore_errors)



    def getUniqueID(self, i_filter_list = []):
        """
        Return a unique ID representing the compile parameters 
        """
        self.makeOptionsConsistent()

        retval = ''

        if not 'compile_misc' in i_filter_list:

            if not 'compile_plane' in i_filter_list:
                if self.plane_spectral_space == 'enable':
                    retval+='_plspec'

                    if self.plane_spectral_dealiasing == 'enable':
                        retval+='_pldeal'

            if not 'compile_sphere' in i_filter_list:
                if self.sphere_spectral_space == 'enable':
                    retval+='_spspec'

                    if self.sphere_spectral_dealiasing == 'enable':
                        retval+='_spdeal'

            if self.gui == 'enable':
                retval+='_gui'

            if self.quadmath == 'enable':
                retval+='_quadmath'

            if self.libfft == 'enable':
                retval+='_fft'


        if self.benchmark_timings == 'enable':
            retval+='_benchtime'

        if not 'compile.parallelization' in i_filter_list:
            if self.sweet_mpi == 'enable':
                retval+='_mpi'

            if self.threading in ['omp']:
                retval+='_th'+self.threading
                
            if self.rexi_thread_parallel_sum == 'enable':
                retval+='_rxthpar'

            if self.rexi_timings_additional_barriers == 'enable':
                retval+='_rxtbar'

            if self.rexi_allreduce == 'enable':
                retval+='_redall'

        retval += '_'+self.mode

        if retval != '':
            retval = 'COMP'+retval

        return retval




    def getOptionList(self):
        return self.__dict__




if __name__ == "__main__":

    p = JobCompileOptions()

    s = p.getSConsParams()
    p.info(s)

    s = p.getProgramName(True)
    p.info(s)

    p.info("FIN")
