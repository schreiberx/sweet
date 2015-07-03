#! /usr/bin/python

import os

import commands
import re
import sys
import platform
import subprocess

from python_mods import CompileXMLOptions

hostname=subprocess.check_output('hostname')
hostname = hostname.replace("\n", '')
hostname = hostname.replace("\r", '')

env = Environment()


###################################################################
# fix environment vars (not imported by default)
###################################################################

# import all environment variables (expecially PATH, LD_LIBRARY_PATH, LD_PATH)
env.Append(ENV=os.environ)


files = os.listdir('src/examples/')
files = sorted(files)

example_programs = []
for f in files:
#	tag = 'example_'
#	if f[0:len(tag)] == tag:
	example_programs.append(f[0:-4])

env['example_programs'] = example_programs


files = os.listdir('src/unit_tests/')
files = sorted(files)

unit_tests_programs = []
for f in files:
	unit_tests_programs.append(f[0:-4])

env['unit_tests_programs'] = unit_tests_programs



################################################################
# XML CONFIGURATION
################################################################


AddOption(	'--xml-config',		dest='xml_config',	type='string',	default='',	help='xml configuration file with compile options')
env['xml_config'] = GetOption('xml_config')

if env['xml_config'] != '':
	env = CompileXMLOptions.load(env['xml_config'], env)


################################################################
# Compile options
################################################################


AddOption(      '--mode',
		dest='mode',
		type='choice',
		choices=['debug', 'release'],
		default='release',
		help='specify compiler to use: debug, release [default: %default]'
)
env['mode'] = GetOption('mode')


AddOption(	'--compiler',
		dest='compiler',
		type='choice',
		choices=['gnu', 'intel', 'llvm', 'pgi'],
		default='gnu',
		help='specify compiler to use: gnu, intel, llvm, pgi [default: %default]'
)
env['compiler'] = GetOption('compiler')


d = ''

if hostname[0:4] == 'mac-':
	if env['compiler'] == 'intel':
		d = ''
	elif env['compiler'] == 'llvm':
		d = '/lrz/sys/compilers/gcc/4.7.3/'

elif hostname == 'martinium':
	if env['compiler'] == 'intel':
		d = 'g++-4.6'

AddOption(	'--gxx-toolchain',
		dest='gxx-toolchain',
		type='string',
		default=d,
		help='specify gcc toolchain for intel and llvm compiler, e.g. g++-4.6, default: deactivated'
)
env['gxx_toolchain'] = GetOption('gxx-toolchain')



AddOption(	'--spectral-space',
		dest='spectral_space',
		type='choice',
		choices=['enable', 'disable'],
		default='disable',
		help="Use spectral space for operations such as folding [default: %default]"
)
env['spectral_space'] = GetOption('spectral_space')

#
# LIB XML
#
AddOption(	'--libxml',
		dest='libxml',
		type='choice',
		choices=['enable','disable'],
		default='disable',
		help='Compile with libXML related code: enable, disable [default: %default]'
)
env['libxml'] = GetOption('libxml')

env.Append(CXXFLAGS=' -DCONFIG_ENABLE_LIBXML='+('1' if env['libxml']=='enable' else '0'))
if env['libxml'] == 'enable':
	env.ParseConfig("xml2-config --cflags --libs")



AddOption(	'--gui',
		dest='gui',
		type='choice',
		choices=['enable','disable'],
		default='disable',
		help='gui: enable, disable [default: %default]'
)
env['gui'] = GetOption('gui')

AddOption(      '--compile-program',
		dest='compile-program',
		type='choice',
		choices=example_programs+unit_tests_programs,
		default='',
		help='Specify program to compile: '+', '.join(example_programs+unit_tests_programs)+' '*80+' [default: %default]'
)
env['compile_program'] = GetOption('compile-program')



env['fortran_source'] = 'disable'


###########################################
# Compile options
###########################################

# set default to first example program
if env['compile_program'] == '':
	env['compile_program'] = example_programs[0]

#
# PROGRAM's NAME
#
if env['compile_program'] in env['unit_tests_programs']:
	env['program_name'] = 'unit_test_'+env['compile_program']

elif env['compile_program'] in env['example_programs']:
	env['program_name'] = 'example_'+env['compile_program']

else:
	print 'FATAL ERROR'
	sys.exit(-1)


env.Append(CXXFLAGS = ' -DSWEET_PROGRAM_NAME='+env['program_name'])



exec_name=env['program_name']

exec_name+='_spectral'+str(env['spectral_space'])

if env['spectral_space'] == 'enable':
	env.Append(CXXFLAGS = ' -DSWEET_USE_SPECTRAL_SPACE=1')
	env.Append(LINKFLAGS=' -lfftw3 -lfftw3_omp')
else:
	env.Append(CXXFLAGS = ' -DSWEET_USE_SPECTRAL_SPACE=0')

env.Append(LINKFLAGS=' -fopenmp')
env.Append(CXXFLAGS=' -fopenmp')

#if env['program_binary_name_override'] != '':
#	exec_name = env['program_binary_name_override']



if env['compiler'] == 'gnu':
	reqversion = [4,6,1]

	#
	# get gcc version using -v instead of -dumpversion since SUSE gnu compiler
	# returns only 2 instead of 3 digits with -dumpversion
	#
	gccv = commands.getoutput('g++ -v').splitlines()

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
		search_string = 'gcc-Version '
		for l in gccv:
			if l[:len(search_string)] == search_string:
				found_line = l
				found = True
	
				gccversion = found_line.split(' ')[1].split('.')
				break

	if not found:
		print search_string+" not found"
		sys.exit(-1)

	for i in range(0, 3):
		if (int(gccversion[i]) > int(reqversion[i])):
			break
		if (int(gccversion[i]) < int(reqversion[i])):
			print 'At least GCC Version 4.6.1 necessary.'
			Exit(1)

	# eclipse specific flag
	env.Append(CXXFLAGS=' -fmessage-length=0')

	# c++0x flag
	env.Append(CXXFLAGS=' -std=c++0x')

	# be pedantic to avoid stupid programming errors
#	env.Append(CXXFLAGS=' -pedantic')

	# SSE 4.2
	env.Append(CXXFLAGS=' -msse4.2')

	# speedup compilation - remove this when compiler slows down or segfaults by running out of memory
	env.Append(CXXFLAGS=' -pipe')

	# activate gnu C++ compiler

	if env['fortran_source']=='enable':
		env.Replace(FORTRAN='gfortran')
		env.Replace(F90='gfortran')
		env.Append(FORTRANFLAGS=' -cpp')
		env.Append(F90FLAGS=' -cpp')
		env.Append(LIBS=['gfortran'])



elif env['compiler'] == 'intel':
	reqversion = [12,1]
	iccversion_line = commands.getoutput('icpc -dumpversion')

	if iccversion_line != 'Mainline':
		iccversion = iccversion_line.split('.')
		for i in range(0, 2):
			if (int(iccversion[i]) > int(reqversion[i])):
				break
			if (int(iccversion[i]) < int(reqversion[i])):
				print 'ICPC Version 12.1 necessary.'
				Exit(1)

	# override g++ called by intel compiler to determine location of header files
	if env['gxx_toolchain'] != '':
		env.Append(CXXFLAGS=' -gxx-name='+env['gxx_toolchain'])
		env.Append(LINKFLAGS=' -gxx-name='+env['gxx_toolchain'])

	env.Append(LINKFLAGS=' -shared-intel')
	env.Append(LINKFLAGS=' -shared-libgcc')
	env.Append(LINKFLAGS=' -debug inline-debug-info')

	# eclipse specific flag
	env.Append(CXXFLAGS=' -fmessage-length=0')

	# c++0x flag
	env.Append(CXXFLAGS=' -std=c++0x')

	# output more warnings
	env.Append(CXXFLAGS=' -w1')

	# compiler option which has to be appended for icpc 12.1 without update1
#	env.Append(CXXFLAGS=' -U__GXX_EXPERIMENTAL_CXX0X__')


	# UBUNTU FIX for i386 systems
	lines = commands.getoutput('uname -i').splitlines()

	for i in lines:
		if i == 'i386':
			env.Append(CXXFLAGS=' -I/usr/include/i386-linux-gnu/')

	# SSE 4.2
	env.Replace(CXX = 'icpc')
	env.Append(CXXFLAGS=' -msse4.2')
	env.Replace(LINK='icpc')

	if env['fortran_source'] == 'enable':
		env.Append(LINKFLAGS=' -lifcore')
		env.Replace(FORTRAN='ifort')
		env.Replace(F90='ifort')
		env.Append(FORTRANFLAGS=' -fpp')
		env.Append(F90FLAGS=' -fpp')


elif env['compiler'] == 'pgi':
	# activate pgi
	env.Replace(CXX = 'pgc++')

	# c++0x flag
	env.Append(CXXFLAGS=' --c++11')

	# last gcc compatible version: 4.6
#	gnu_compat_version = commands.getoutput("`g++-4.6 -v 2>&1 | grep ' version ' | cut -d' ' -f3")

	# set gnu version
#	env.Append(CXXFLAGS=' --gnu_version '+gnu_compat_version)
#	env.Append(CXXFLAGS=' --gnu_version 4.8.0')
	env.Append(CXXFLAGS=' -D__GXX_EXPERIMENTAL_CXX0X__ --gnu_extensions')

	# be pedantic to avoid stupid programming errors
#	env.Append(CXXFLAGS=' -pedantic')

	# SSE 4.2
#	env.Append(CXXFLAGS=' -msse4.2')

	# speedup compilation - remove this when compiler slows down or segfaults by running out of memory
#	env.Append(CXXFLAGS=' -pipe')

	if env['fortran_source'] == 'enable':
		print "TODO: PGI compiler not yet supported with fortran enabled"
		Exit(-1)


elif env['compiler'] == 'llvm':
	reqversion = [3,1]

	if env['gxx_toolchain'] != '':
		env.Append(CXXFLAGS=' --gcc-toolchain='+env['gxx_toolchain'])

	# eclipse specific flag
	env.Append(CXXFLAGS=' -fmessage-length=0')

	# c++0x flag
	env.Append(CXXFLAGS=' -std=c++0x')

	# be pedantic to avoid stupid programming errors
	env.Append(CXXFLAGS=' -pedantic')

	# SSE 4.2
	env.Append(CXXFLAGS=' -msse4.2')

	# speedup compilation - remove this when compiler slows down or segfaults by running out of memory
	env.Append(CXXFLAGS=' -pipe')

	# activate gnu C++ compiler

	# todo: fix me also for intel mpicxx compiler
	env.Replace(CXX = 'clang++')

	if env['fortran_source'] == 'enable':
		print "TODO: LLVM compiler not yet supported with fortran enabled"
		Exit(-1)



if env['mode'] == 'debug':
	env.Append(CXXFLAGS=' -DSWEET_DEBUG_MODE=1')

	if env['compiler'] == 'gnu':
		env.Append(CXXFLAGS=' -O0 -g3 -Wall')

		# integer overflow check
		env.Append(CXXFLAGS=' -ftrapv')

	elif env['compiler'] == 'llvm':
		env.Append(CXXFLAGS=' -O0 -g3 -Wall')

	elif env['compiler'] == 'intel':
		env.Append(CXXFLAGS=' -O0 -g -traceback')
#		env.Append(CXXFLAGS=' -fp-trap=common')

	elif env['compiler'] == 'pgi':
		env.Append(CXXFLAGS=' -O0 -g -traceback')


	if env['fortran_source'] == 'enable':
		if env['compiler'] == 'gnu':
			env.Append(FORTRANFLAGS=' -g -O0')
			env.Append(F90FLAGS=' -g -O0')
		elif env['compiler'] == 'intel':
			env.Append(FORTRANFLAGS=' -g -O0 -traceback')
			env.Append(F90FLAGS=' -g -O0 -traceback')


elif env['mode'] == 'release':
	env.Append(CXXFLAGS=' -DSWEET_DEBUG_MODE=0')

	if env['compiler'] == 'gnu':
		env.Append(CXXFLAGS=' -O3 -mtune=native')

	if env['compiler'] == 'llvm':
		env.Append(CXXFLAGS=' -O4 -mtune=native')

	elif env['compiler'] == 'intel':
		env.Append(CXXFLAGS=' -O3 -fno-alias')
		env.Append(CXXFLAGS=' -ipo')
		env.Append(CXXFLAGS=' -fast')

	elif env['compiler'] == 'pgi':
		env.Append(CXXFLAGS='-O3 -fast -Mipa=fast,inline -Msmartalloc')

	if env['fortran_source'] == 'enable':
		if env['compiler'] == 'gnu':
			env.Append(FORTRANFLAGS=' -O2')
			env.Append(F90FLAGS=' -O2')
		elif env['compiler'] == 'intel':
			env.Append(FORTRANFLAGS=' -O3')
			env.Append(F90FLAGS=' -O3')
#			env.Append(FORTRANFLAGS=' -ipo')
#			env.Append(F90FLAGS=' -ipo')
#			env.Append(FORTRANFLAGS='-fast')
#			env.Append(F90FLAGS=' -fast')




if env['gui'] == 'enable':
	# compile flags
	env.Append(CXXFLAGS=' -I'+os.environ['HOME']+'/local/include')
	env.Append(CXXFLAGS=' -DSWEET_GUI=1')

	# linker flags

	# add nvidia lib path when running on atsccs* workstation
#	hostname = commands.getoutput('uname -n')
#	if re.match("atsccs.*", hostname) or  re.match("laptop.*", hostname):
#		env.Append(LIBPATH=['/usr/lib/nvidia-current/'])

#	env.Append(LIBPATH=[os.environ['HOME']+'/local/lib'])
	env.Append(LIBS=['GL'])

	reqversion = [2,0,0]
	sdlversion = commands.getoutput('sdl2-config --version').split('.')

	for i in range(0, 3):
		if (int(sdlversion[i]) > int(reqversion[i])):
			break;
		if (int(sdlversion[i]) < int(reqversion[i])):
			print 'libSDL Version 2.0.0 necessary.'
			Exit(1)

	env.ParseConfig("sdl2-config --cflags --libs")
	env.ParseConfig("pkg-config SDL2_image --cflags --libs")
	env.ParseConfig("pkg-config freetype2 --cflags --libs")
else:
	env.Append(CXXFLAGS=' -DSWEET_GUI=0')



#
# build directory
#
#build_dir='build/build_'+exec_name+'/'
build_dir='/tmp/scons_build_'+exec_name+'/'


env.Append(CPPPATH = ['/usr/local/include', '/usr/include'])
# also include the 'src' directory to search for dependencies
env.Append(CPPPATH = ['.', 'src/'])


######################
# get source code files
######################


env.src_files = []

Export('env')
SConscript('src/SConscript', variant_dir=build_dir, duplicate=0)
Import('env')

print
print '            Program: '+env['program_name']
print 'Building executable: '+exec_name
print

env.Program('build/'+exec_name, env.src_files)

