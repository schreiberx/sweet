#! /usr/bin/env python3

from SWEET import *
from mule.postprocessing.JobData import *
from mule.postprocessing.JobsDataConsolidate import *
import glob
import os
import re
import sys



class JobsData(InfoError):

	def __init__(
		self,
		job_dirs = 'job_bench_*',
		verbosity = 10,
	):
		"""
		Load the data of all jobs in the current folder

		Parameters:
		-----------
			job_dirs: string / list
				if string:
					Pattern with wildcards to autodetect job directories, e.g. 'job_bench_*'
		"""
		InfoError.__init__(self, 'JobsData')

		self.verbosity = verbosity

		# If job_dirs is just a string, it's a search pattern
		if isinstance(job_dirs, str):
			job_dirs = glob.glob(job_dirs)
		elif isinstance(job_dirs, list):
			pass
		else:
			raise Exception("Unknown type for job_dirs")

		# Load raw job data
		self.__load_job_raw_data(job_dirs)

		# Create flattened data
		self.__create_flattened_data()




	def __load_job_raw_data(
			self,
			job_dirs = []
	):
		"""
		Parse all output.out files and extract all kind of job output information

		Return a dictionary with content from the job directories
		{
			[name of job directory] :
			{
				#
				# Dictionary with data from job generation
				# (read from [jobdir]/jobgeneration.pickle)
				#
				'jobgeneration':
				{
					'compile': [...],
					'runtime': [...],
					'parallelization': [...],
					'platforms_platform': [...],
					'platform_resources': [...],
				},
				'output':
				{
					'SimulationBenchmarkTimings.main': [value],
					'SimulationBenchmarkTimings.main_simulationLoop': [value],
					[...]
				}
			}
		"""

		self.__jobs_data = {}
		for job_dir in job_dirs:
			if self.verbosity > 5:
				self.info("")
				self.info("Processing '"+job_dir+"'")


			job = JobData(job_dir, verbosity=self.verbosity)
			if 'jobgeneration.job_unique_id' not in job.get_flattened_data():
				# Be backward compatible
				self.__jobs_data[job_dir] = job
			else:
				job_unique_id = job.get_flattened_data()['jobgeneration.job_unique_id']
				self.__jobs_data[job_unique_id] = job



	def get_jobs_data(self):
		"""
		Return the raw job information data
		Warning: This storage format is likely to change!
		"""
		return self.__jobs_data



	def __create_flattened_data(self):
		"""
		Return a dictionary with a flattened hierarchy

		This joins all hierarchical structures with a '.'
		
		E.g. the value
			self.jobs_raw_data['jobgeneration'].parallelization.num_cores_per_socket
		is stored in
			retval['jobgeneration.parallelization.num_cores_per_socket']
		"""

		self.__flattened_jobs_data = {}

		# For all jobs
		for job_key, job_data in self.__jobs_data.items():
			self.__flattened_jobs_data[job_key] = job_data.get_flattened_data()



	def get_flattened_data(self):
		return self.__flattened_jobs_data



if __name__ == "__main__":

	verbosity = 10
	verbosity = 0

	if True:
		j = JobsData(verbosity=verbosity)
#		c = JobsDataConsolidate(verbosity=verbosity)
		d = j.get_flattened_data()

		for jobdir, job_data in d.items():
			print("*"*80)
			print("Data for '"+jobdir+"'")
			print("*"*80)
			for key, value in job_data.items():
				print(key+" => "+str(value))

	sys.exit(1)

	if False:
		j = JobsData(verbosity=verbosity)

		jobs_raw_data = j.get_jobs_raw_data()
		for key, values in jobs_raw_data.items():
			print(key)
			print(values)
			print("")

	#if True:
	if False:
		j = JobsData(verbosity=verbosity)
		d = j.get_flattened_data()

		for jobdir, job_data in d.items():
			print("")
			print("Data for '"+jobdir+"'")
			for key, value in job_data.items():
				print(key+" => "+str(value))

	#if True:
	if False:
		j = JobsData(verbosity=verbosity)
		d = j.get_flattened_data()

		# Just print the first item
		#if False:
		if True:
			for jobdir, job_data in d.items():
				print("")
				print("Data for '"+jobdir+"'")
				for key, value in job_data.items():
					print(key+" => "+str(value))
				break

		g = j.create_groups(['runtime.timestepping_method'])

		#for group_id, group_jobs in g.items():
		for group_id in sorted(groups):
			group_jobs = groups[group_id]

			print(group_id)
			for jobdir, jobdata in group_jobs.items():
				print(" + "+jobdir+"\t"+jobdata['output.benchmark_timings.main_simulationloop'])


	#if True:
	if False:
		j = JobsData(verbosity=verbosity)

		data_table = j.create_data_table(
				['runtime.timestepping_method'],
				'parallelization.num_threads_per_rank',
				'output.benchmark_timings.main_simulationloop'
			)

		print("Data table:")
		j.print_data_table(data_table)

