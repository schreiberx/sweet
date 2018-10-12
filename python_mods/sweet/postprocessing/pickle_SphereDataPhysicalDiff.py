#! /usr/bin/env python3

import sys
import math
import glob

from .JobData import *
from .JobsData import *
from .SphereDataPhysicalDiff import *



class pickle_SphereDataPhysicalDiff:

	def __init__(
			self,
			ref_file_ending,
			jobdir_pattern = None,
		):
		"""
		Generate the .pickle files in each job directory based on simulation output files given by 'ref_file_ending'

		Parameters
		----------
		ref_file_ending: str
			string with ending of reference files
		"""
		self._setup(ref_file_ending, jobdir_pattern)
		pass



	def get_matching_files_in_dir(
		self,
		path,
		needle
	):
		files = os.listdir(path)
		retval = []
		for f in files:
			if needle in f:
				retval.append(f)
		return retval


	def _setup(
			self,
			ref_file_ending = None,
			jobdir_pattern = None,
		):

		if jobdir_pattern == None:
			jobdir_pattern = './job_bench*'

		if ref_file_ending == None:
			ref_file_ending = "_t00000000120.00000000.csv"

		j = JobsData(jobdir_pattern, verbosity=0)
		jobs = j.get_flattened_data()

		no_reference_job_unique_id_found = True

		if len(jobs) == 0:
			raise Exception("No jobs found!")

		for key, job in jobs.items():
			print("Processing "+key)

			# Sort out jobs which don't have a reference job id
			# These jobs are likely the reference jobs themselves
			if 'jobgeneration.reference_job_unique_id' not in job:
				continue

			no_reference_job_unique_id_found = False

			reference_job_unique_id = job['jobgeneration.reference_job_unique_id']
			print(" + ref job id: "+reference_job_unique_id)

			ref_key = None
			for skey, sjob in jobs.items():
				if sjob['jobgeneration.job_unique_id'] == reference_job_unique_id:
					ref_key = skey

			if ref_key == None:
				print("Fatal: missing reference job with id "+reference_job_unique_id)
				raise Exception("Reference job not found!")

			# Load reference job
			ref_job = jobs[ref_key]

			# Load reference files
			ref_files = self.get_matching_files_in_dir(ref_job['jobgeneration.job_dirpath'], ref_file_ending)

			for ref_file in ref_files:
				s = None
				try:
					s = SphereDataPhysicalDiff(
							ref_job['jobgeneration.job_dirpath']+'/'+ref_file,
							job['jobgeneration.job_dirpath']+'/'+ref_file
					)
				except Exception as e:
					print("Error occured which is ignored")
					print(str(e))
					# Ignore missing files
					continue

				s.print()

				pickle_filename = 'sphere_data_diff_'+ref_file.replace('output_', '').replace(ref_file_ending, '')+'.pickle'
				print("Writing file "+pickle_filename)
				s.write_file(job['jobgeneration.job_dirpath']+'/'+pickle_filename)

			print(ref_key)
			print("")



		if no_reference_job_unique_id_found:
			print("*"*80)
			print("Warning: No data generated")
			print("No job with a reference_job_unique_id found!")
			print("*"*80)



if __name__ == '__main__':

	p = pickle_SphereDataPhysicalDiff()

