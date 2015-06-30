#
#  Created on: Dec 7, 2014
#      Author: Martin Schreiber <schreiberx@gmail.com>
#
# Unit test support module
#

import subprocess
import time
import os
import re
import sys

#
# This class abstracts the unit tests
#
class UnitTests:

	# name of main unit test
	main_unit_test = None

	# subtest id
	subtest_id = None
	subtest_id_max = None

	# previous compiler command to avoid duplicated compilation
	# e.g. if only different simulation parameters are used
	prev_compiler_cmd=''

	# how many processes to use for compilation
	parallel_compilation_processes = 4


	# constructor with argv program parameters
	def __init__(self, argv):
		# 1st parameter: main unit test
		self.main_unit_test = None
		if len(argv) > 1:
			self.main_unit_test = argv[1]

		self.subtest_id = None
		self.subtest_max = None
		if len(argv) > 2:
			self.subtest_id = int(argv[2])
		if len(argv) > 3:
			self.subtest_id_max = int(argv[3])


	def testCompiler(
			self,
			command_list,
			reqversion,
		):
		print 'Testing compiler '+command_list[0]+' ... ',

		try:
			output = subprocess.check_output(command_list)
		except OSError as e:
			print 'Failed to execute command "'+command_list[0]+'":'
			print e
			sys.exit(-1)
		except subprocess.CalledProcessError as e:
			print 'Failed to execute '+command_list[0]+':'
			print e
			sys.exit(-1)

		lines = output.splitlines()
		found = False

		# extract line with program version
		version_line = lines[0]
		if len(version_line) == 0:
			version_line = lines[1]

		# subversion numbers found
		numbers_hit = -1
	
		# iterate over version numbers consisting out of 4, 3, and 2 numbers to find the longest one
		for numbers in range(4, 1, -1):
			# support formats such as
			# 1.2.3 (intel, gnu, clang)
			# 1.2-3 (pgi)
			matchstring = "(\d+)"
			for i in range(1, numbers):
				matchstring = matchstring+"[.-](\d+)"

			m = re.search(matchstring, version_line)

			if m != None:
				numbers_hit = numbers
				break

		if numbers_hit == -1:
			print 'Version number not found!'
			sys.exit(-1)

		#print 'Matching found with '+str(numbers_hit)+' numbers of subversions and matching string '+matchstring

		version_subs = m.groups()
		#print version_subs

		for i in range(0, min(len(reqversion), numbers_hit)):
			if (int(version_subs[i]) > reqversion[i]):
				break
			if (int(version_subs[i]) < reqversion[i]):
				print 'ERROR: Requested version '+('.'.join(map(str, reqversion)))+' for '+command_list[0]+' not found'
				sys.exit(-1)

		print 'OK'


	def testCompilers(self, compiler_list):
		for c in compiler_list:
			if c == 'gnu':
				self.testCompiler(['g++', '--version'], [4,6,1])

			elif c == 'intel':
				self.testCompiler(['icpc', '--version'], [12,1])

			elif c == 'pgi':
				self.testCompiler(['pgc++', '--version'], [14,10,0])

			elif c == 'llvm':
				self.testCompiler(['clang++', '--version'], [3,5,0])
				


	#
	# executed the unit tests given by i
	#
	def unitTest(
		self,
		name,		# name of main test
		subtest_id,	# id of subtest
		i,		# tuple [compiler command, execution command, expected return value [default: 0]]
		outputAll = False
	):
		expected_returncode = i[2] if len(i) == 3 else 0
		logfile = '/tmp/tests_'+name+'_'+str(subtest_id)+'.log'

		cmd_append=' 2>&1 >'+logfile

		CRED = '\033[91m'
		CGREEN = '\033[92m'
		CDEFAULT = '\033[0m'

		def print_err(s):
			print CRED+s+CDEFAULT

		def print_ok(s):
			print CGREEN+s+CDEFAULT

		print "Running "+str(name)+' / '+str(subtest_id)

		compiler_cmd = i[0]+' -j '+str(self.parallel_compilation_processes)

		if compiler_cmd != self.prev_compiler_cmd:
			print "COMPILING: "+compiler_cmd
			p = subprocess.Popen(['make clean'+cmd_append], shell=True)
			p.wait()
			if p.returncode != 0:
				print_err(" > FAILED TO MAKE CLEAN!")
				return

			self.prev_compiler_cmd = compiler_cmd

			p = subprocess.Popen([compiler_cmd+cmd_append], shell=True)
			p.wait()
			if p.returncode != 0:
				print_err(" > FAILED TO COMPILE! ("+compiler_cmd+")")
				sys.exit(-1)
		else:
			print "Using previously compiled program"

		startTime = time.time()

		build_dir='./build/'
		files = os.listdir(build_dir)
		match=None
		for j in files:
			print j
			if not os.path.isdir(build_dir+j):
#				print "not isdir "+j
				if match != None:
					print "Two executables found, conflicting builds"
					sys.exit(-1)

				match=build_dir+j

		if match==None:
			print "No executable found!"
			sys.exit(-1)

		exec_cmd = match
		print "FOUND EXECUTABLE: "+exec_cmd

		exec_cmd_and_params=(i[1]%exec_cmd)+' '+cmd_append

		print "EXECUTING: "+exec_cmd_and_params

		p = subprocess.Popen([exec_cmd_and_params], shell=True)
		p.wait()

		if p.returncode != expected_returncode:
			print_err(" > TEST FAILED with return code "+str(p.returncode))
			print_err("See "+logfile+" for errors")
			sys.exit(-1)
		else:
			print_ok(" > TEST OK ("+str(time.time() - startTime)+" Seconds)")

			cps_str = 'Cells per second (CPS): '
			with open(logfile) as f:
				for l in f:
					if outputAll:
						print l,
					else:
						if l[0:len(cps_str)] == cps_str:
							cps = l[len(cps_str):]
							print_ok(exec_cmd+"\tCPS="+cps)
							break

		sys.stdout.flush()
		sys.stderr.flush()

		print "*"*60



	def run(self, name, tests, outputAll = False):

		if self.main_unit_test != None:
			if name != self.main_unit_test:
				return

		print
		print "*"*60
		print "* MAIN TEST NAME: "+name
		print "*"*60

		if self.subtest_id_max != None:
			for i in range(self.subtest_id, min(len(tests), self.subtest_id_max)):
				self.unitTest(name, i, tests[i], outputAll)
			sys.exit(0)

		if self.subtest_id != None:
			self.unitTest(name, self.subtest_id, tests[self.subtest_id], outputAll)
			sys.exit(0)

		for i in range(0, len(tests)):
			self.unitTest(name, i, tests[i], outputAll)

