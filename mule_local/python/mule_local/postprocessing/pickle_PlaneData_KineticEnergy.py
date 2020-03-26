#! /usr/bin/env python3

import sys
import math
import glob

from mule.InfoError import *
from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from mule_local.postprocessing.PlaneData_KineticEnergy import *



class pickle_PlaneData_KineticEnergy(InfoError):

	def __init__(
			self,
			ref_file_ending = None,
			jobdir_pattern = None,
			job_dirs = None,
			params = [],
		):
		"""
		Generate the .pickle files in each job directory based on simulation output files given by 'ref_file_ending'

		Parameters
		----------
		ref_file_ending: str
			string with ending of reference files

		jobdir_pattern: str
			string with matching pattern for job directories

		job_dirs: list
			list with job directories

		params: list of strings
			list with optional parameters

			'ignore_missing_file':
				Don't throw an error if a file is missing

			'only_last_file':
				Process only last file

			ohter parameters:
				see PlaneData_KineticEnergy.py
		"""

		InfoError.__init__(self, "pickle_PlaneData_KineticEnergy")
		self._setup(ref_file_ending=ref_file_ending, jobdir_pattern=jobdir_pattern, job_dirs=job_dirs, params=params)



	def _setup(
			self,
			ref_file_ending = None,
			jobdir_pattern = None,
			job_dirs = None,
			params = [],
		):

		if job_dirs != None:
			j = JobsData(job_dirs=job_dirs, verbosity=0)

		else:
			if jobdir_pattern == None:
				jobdir_pattern = './job_bench*'

			j = JobsData(jobdir_pattern, verbosity=0)

		jobs = j.get_flattened_data()

		no_reference_job_unique_id_found = True

		if len(jobs) == 0:
			raise Exception("No jobs found!")

		for key, job in jobs.items():
			print("Processing "+key)

			# job directory path
			job_dirpath = job['jobgeneration.job_dirpath']

			# u-velocity 
			u_vel_files = glob.glob(job_dirpath+'/output_*_u_*.csv')

			if len(u_vel_files) == 0:
				self.info("WARNING")
				self.info("WARNING: No velocity files found")
				self.info("WARNING: However, there should be at least one velocity file (the one at t=0)")
				self.info("WARNING")


			# only process the very last file
			if 'only_last_file' in params:

				u_vel_files.sort()
				u_vel_files = [u_vel_files[-1]]

				# determine (time-depending) ending of reference file
				pos = u_vel_files[0].rfind('_')
				if pos < 0:
					raise Exception("File ending not found for reference file '"+u_vel_files[0]+"'")


			# Iterate over all velocity files for different time stamps
			for u_vel_file in u_vel_files:
				v_vel_file = u_vel_file.replace('_u_', '_v_')

				try:
					s = PlaneData_KineticEnergy(
							u_vel_file,
							v_vel_file,
							params
					)

				except FileNotFoundError as e:
					# Ignoring missing files should be configured via "ignore_missing_file" parameter, see above
					if "ignore_missing_file" in params:
						self.info("Ignoring Error:")
						self.info(str(e))
						continue

					raise Exception(e)

				except IOError as e:
					# Ignoring missing files should be configured via "ignore_missing_file" parameter, see above
					if "ignore_missing_file" in params:
						self.info("Ignoring Error:")
						self.info(str(e))
						continue

					raise Exception(e)

				#s.print()

				if 'only_last_file' in params:
					pickle_filename = 'plane_data_kinetic_energy.pickle'

				else:

					# determine (time-depending) ending of reference file
					pos = u_vel_files[0].rfind('_')
					if pos < 0:
						raise Exception("File ending not found for reference file '"+u_vel_files[0]+"'")

					pickle_file_ending = u_vel_file[pos:]
					pickle_file_ending = pickle_file_ending.replace('.csv', '')
					print("pickle_file_ending: "+pickle_file_ending)

					pickle_filename = 'plane_data_kinetic_energy'+pickle_file_ending+'.pickle'

				print("Writing file "+pickle_filename)
				s.write_file(job['jobgeneration.job_dirpath']+'/'+pickle_filename)



if __name__ == '__main__':

	#p = pickle_PlaneData_KineticEnergy(params=['only_last_file'])
	p = pickle_PlaneData_KineticEnergy()

