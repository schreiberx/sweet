#! /usr/bin/env python3


import os
import sys
import re
import importlib
from inspect import getmembers
from InfoError import *

__all__ = ['SWEETPlatforms']


class SWEETPlatforms(InfoError):

	def __init__(	
			self,
			platform_id_override = None	# Override automatic detection of platform
	):
		"""
		Initialize platform and setup target platform

		Sideeffects:
		-----------
		self.platforms:
			Dictionary with all available platforms
			
		self.platform_id:
			Platform ID which will be used

		self.module:
			Module of platform to be used

		self.functions:
			Dictionary of functions for platform to be used

		self.platform_dirpath:
			Path to directory where platform information can be found		
		"""

		InfoError.__init__(self, "SWEETPlatforms")

		self.sweetroot = os.environ.get('SWEET_ROOT')
		if self.sweetroot == None:
			self.error("SWEET environment variables not loaded!")

		if platform_id_override == None:
			# Check for environment variable SWEET_PLATFORM
			platform_id_override = os.environ.get('SWEET_PLATFORM')

		self.platform_id_override = platform_id_override

		self.platforms_dir = os.path.normpath(self.sweetroot+'/platforms/')

		self.hline()
		self.info("Platforms description path: "+self.platforms_dir)


		# Required interfaces by platforms
		self.interface_names = {
			'get_platform_id',
			'get_platform_autodetect',
			'get_platform_resources',

			# a reference to SWEETJobGeneration is handed over as a parameter to all these functions
			'jobscript_setup',			# setup return of job script content
			'jobscript_get_header',			# header (e.g. scheduler print_information) for job script
			'jobscript_get_exec_prefix',		# prefix before MPI executable
			'jobscript_get_exec_command',		# MPI execution (something like "mpirun -n ### "
			'jobscript_get_exec_suffix',		# suffix aftere MPI executable
			'jobscript_get_footer',			# footer (e.g. postprocessing) for job script

			'jobscript_get_compile_command',	# suffix aftere MPI executable
		}

		self.platform_id = None

		self.p_setup()




	def p_load_module_interfaces(
			self,
			module
	):
		# Import all members from module
		members = getmembers(module)

		# Setup dictionary for found interfaces
		interfaces = {i: None for i in self.interface_names}

		for m in members:
			if m[0] in self.interface_names:
				interfaces[m[0]] = m[1]

		self.info("Loading platform description for '"+ module.__name__+"'")
		for i in interfaces:
			if interfaces[i] == None:
				self.error("Interface '"+i+"' in module '"+module.__name__+"' not found")
		self.info(" + platform_id: '"+ interfaces['get_platform_id']()+"'")

		return interfaces

	"""
	Setup and validation

	self.platforms:	Dictionary with all platforms
	self.platform: Target platform
	"""
	def p_setup(
			self
	):
		self.platforms = {}

		# Iterate over all possible packages
		files = os.listdir(self.platforms_dir)

		# Sort to include priority enumeration, e.g. 50_cheyenne
		files.sort()

		pattern = re.compile('[0-9]{2}_[a-z]+')
		for f in files:
			if f in ['.', '..']:
				continue

			if pattern.match(f) == None:
				continue

			sys.path.append(self.platforms_dir)
			module = importlib.import_module(f+'.SWEETPlatform')
			sys.path.remove(self.platforms_dir)

			interfaces = self.p_load_module_interfaces(module)

			platform_id = interfaces['get_platform_id']()

			if platform_id in self.platforms:
				self.error("Duplicate platform id '"+platform_id+"' detected in modules '"+self.platforms[platform_id].module.__name__+"' and '"+module.__name__+"'")

			self.platforms[platform_id] = {
				'interfaces': interfaces,
				'module': module,
				'dirname': f,
			}

		self.platform = None

		if self.platform_id_override != None:
			if self.platform_id_override not in self.platforms:
				self.error("Platform override error: Platform with id '"+self.platform_id_override+"' not found")

			self.platform = self.platforms[self.platform_id_override]

		else:
			for key, platform in self.platforms.items():
				if platform['interfaces']['get_platform_autodetect']():
					self.platform = platform
					break

		if self.platform == None:
			self.error("No platform available!")


		self.platform_dirpath = os.path.normpath(self.platforms_dir+'/'+self.platform['dirname'])

		self.module = self.platform['module']

		class Functions: pass
		self.functions = Functions
		for i in self.platform['interfaces']:
			setattr(self.functions, i, getattr(self.module, i))

		self.platform_id = self.functions.get_platform_id()

		self.info("Using platform '"+str(self.platform_id)+"'")



if __name__ == "__main__":

	print()
	p = SWEETPlatforms()
	p.info("TEST - Using platform ID: "+p.platform_id)

	print()
	platform_id_override = "coolmuc_mpp2"
	p = SWEETPlatforms(platform_id_override)
	p.info("TEST - Using platform ID: "+p.platform_id)

	if p.platform_id != platform_id_override:
		p.error("TEST - Platform override '"+platform_id_override+"' not working")

	if p.platform_id != p.functions.get_platform_id():
		p.error("TEST - Invalid platform return value of p.functions.get_platform_id()")

