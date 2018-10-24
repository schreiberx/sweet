#! /usr/bin/env python3

import sys
import math
import glob

from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from sweet.postprocessing.PlaneDataPhysicalDiff import *



class pickle_PlaneDataPhysicalDiff:

	def __init__(
			self,
			ref_file_ending = None,
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





	def _setup(
			self,
			ref_file_ending = None,
			jobdir_pattern = None,
		):

		if jobdir_pattern == None:
			jobdir_pattern = './job_bench*'

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
					if ref_key != None:
						raise Exception("Reference job already found and now there's another one? This is probably not what you wanted")

					ref_key = skey

			if ref_key == None:
				print("Fatal: missing reference job with id "+reference_job_unique_id)
				print("Fatal: reference job was intended for job with dirpath: "+job['jobgeneration.job_dirpath'])
				raise Exception("Reference job not found!")

			# Load reference job
			ref_job = jobs[ref_key]

			if 'output.reference_filenames' in sjob:
				# Were the reference filenames provided?
				ref_files = sjob['output.reference_filenames'].split(";")

				# Now we have to find the file ending
				# We guess that this starts at the last '_' character in the filename

				pos = ref_files[0].rfind('_')
				if pos < 0:
					raise Exception("File ending not found for reference file '"+ref_files[0]+"'")

				use_ref_file_ending = ref_files[0][pos:]
				print("use_ref_file_ending: "+use_ref_file_ending)

			else:

				if ref_file_ending != None:
					use_ref_file_ending = ref_file_ending
				else:
					print("*"*80)
					print(ref_job['runtime.simtime'])

					# "output_%s_t%020.8f.csv"
					use_ref_file_ending = "_t{:020.8f}.csv".format(float(ref_job['runtime.simtime'])/(60*60))

				if use_ref_file_ending == "":
					raise Exception("No reference file ending provided / found")

				# Load reference files
				ref_files = []
				files = os.listdir(ref_job['jobgeneration.job_dirpath'])
				for f in files:
					if use_ref_file_ending in f:
						ref_files.append(f)

			if len(ref_files) == 0:
				print("No reference files found!")
				print("*"*80)
				print("Reference directory: "+ref_job['jobgeneration.job_dirpath'])
				print("Reference file endings: "+use_ref_file_ending)
				print("*"*80)
				raise Exception("Reference files not found!")

			for ref_file in ref_files:
				s = None
				try:
					s = PlaneDataPhysicalDiff(
							ref_job['jobgeneration.job_dirpath']+'/'+ref_file,
							job['jobgeneration.job_dirpath']+'/'+ref_file
					)
				except Exception as e:
					print("Error occured which is ignored")
					print(str(e))
					# Ignore missing files
					continue

				s.print()

				pickle_filename = 'plane_data_diff_'+ref_file.replace('output_', '').replace(use_ref_file_ending, '')+'.pickle'
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

	p = pickle_PlaneDataPhysicalDiff()

