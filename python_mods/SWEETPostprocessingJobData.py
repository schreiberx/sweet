#! /usr/bin/env python3

from SWEET import *
import glob
import os
import re
import sys



class SWEETPostprocessingJobData(InfoError):

	def __init__(
		self,
		jobdir = None,
		verbosity = 10
	):
		InfoError.__init__(self, 'SWEETPostprocessingJobData')

		self.verbosity = verbosity

		# Load raw job data
		self.__load_job_raw_data(jobdir)

		# Create flattened data
		self.__create_flattened_data()



	def __parse_job_output(
		self,
		output_lines
	):
		retdict = {}
		for l in output_lines:
			m=re.match("^\[MULE\] ([^ :]*): (.*)$", l)
			if m != None:
				tag = m.group(1)
				data = m.group(2)
				retdict[tag] = data

		return retdict



	def __load_job_raw_data(
			self,
			jobdir = None
	):
		"""
		Parse all output.out files and extract all kind of job output information

		Return a dictionary with content from the job directories
		{
			#
			# Dictionary with data from job generation
			# (read from [jobdir]/jobgeneration.pickle)
			#
			'jobgeneration':	# From jobgeneration.pickle
			{
				'compile': [...],
				'runtime': [...],
				'parallelization': [...],
				'platforms_platform': [...],
				'platform_resources': [...],
			},
			'output':		# From output.out with prefix [MULE]
			{
				'simulation_benchmark_timings.main': [value],
				'simulation_benchmark_timings.main_simulationLoop': [value],
				[...]
			},
			'[filename(.pickle)]':
			{
				[...]
			}
			},
			'[filename(.pickle)]':
			{
				[...]
			}
		"""

		# 
		self.__job_raw_data = {}

		if self.verbosity > 5:
			self.info("")
			self.info("Processing '"+jobdir+"'")

		"""
		* Process 'output.out'
		"""
		if self.verbosity > 5:
			self.info("Loading job output file 'output.out'")
		outfile = jobdir+'/output.out'
		#try:
		if True:
			with open(outfile, 'r') as f:
				content = f.readlines()
				self.__job_raw_data['output'] = self.__parse_job_output(content)
		"""
		except Exception as err:
			print("*"*80)
			print("* ERROR opening '"+outfile+"' (ignoring)")
			print("* "+str(err))
			print("*"*80)
			continue

"""


		"""
		Process 'jobgeneration.pickle'
		"""
		# Ensure that 'jobgeneration.pickle' exists
		jobgenerationfile = jobdir+'/jobgeneration.pickle'
		if self.verbosity > 5:
			self.info("Loading 'jobgeneration.pickle'")
		j = SWEETJobGeneration(dummy_init=True)
		self.__job_raw_data['jobgeneration'] = j.load_attributes_dict(jobgenerationfile)

		"""
		Process other '*.pickle'
		"""
		pickle_files = glob.glob(jobdir+'/*.pickle')

		# Iterate over all found pickle files
		for picklefile in pickle_files:
			filename = os.path.basename(picklefile)
			tag = filename.replace('.pickle', '')
			if tag == 'jobgeneration':
				continue

			if self.verbosity > 5:
				self.info("Loading pickle file '"+filename+"'")

			import pickle
			with open(picklefile, 'rb') as f:
				self.__job_raw_data[tag] = pickle.load(f)



	def get_job_raw_data(self):
		"""
		Return the raw job information data
		Warning: This storage format is likely to change!
		"""
		return self.__job_raw_data



	def __create_flattened_data(self):
		"""
		Return a dictionary with a flattened hierarchy

		This joins all hierarchical structures with a '.'
		
		E.g. the value
			self.__job_raw_data['jobgeneration'].parallelization.num_cores_per_socket
		is stored in
			self.__flattened_job_data['jobgeneration.parallelization.num_cores_per_socket']
		"""

		self.__flattened_job_data = {}

		# For all raw data
		for raw_key, raw_value in self.__job_raw_data.items():
			if raw_key == 'output':
				for key, data in raw_value.items():
					key = 'output.'+key.lower()
					self.__flattened_job_data[key] = data

			elif raw_key == 'jobgeneration':
				for group_name, group_value in raw_value.items():
					if isinstance(group_value, dict):
						for key, attr in group_value.items():
							key = group_name+"."+key
							self.__flattened_job_data[key] = attr
					elif isinstance(group_value, (int, float, str)):
						self.__flattened_job_data[group_name] = group_value
					else:
						raise Exception("Unknown type in pickled data")

			else:
				for key, data in raw_value.items():
					key = key.lower()
					self.__flattened_job_data[raw_key+'.'+key] = data
	


	def get_flattened_data(self):
		return self.__flattened_job_data



if __name__ == "__main__":

	verbosity = 0

	jobdirs = []

	if len(sys.argv) > 1:
		jobdirs = sys.argv[1:]
	else:
		print("")
		print("Usage:")
		print("	"+sys.argv[0]+" [jobdir 1] [jobdir 2] ...")
		print("")
		sys.exit(1)

	for j in jobdirs:

		if verbosity > 5:
			print("*"*80)
			print("Job directory: "+j)
			print("*"*80)

		j = SWEETPostprocessingJobData(jobdir = j, verbosity=verbosity)
		d = j.get_flattened_data()

		for key, value in d.items():
			print(key+" => "+str(value))
	#			break
