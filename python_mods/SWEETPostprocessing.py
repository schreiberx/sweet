#! /usr/bin/env python3

from SWEET import *
import glob
import os
import re
import sys



class SWEETPostprocessing(InfoError):

	def __init__(
		self,
		jobs_dir = None,
		verbosity = 10
	):
		InfoError.__init__(self, 'SWEETPostprocessing')

		self.verbosity = verbosity

		# Load raw job data
		self.__load_job_raw_data(jobs_dir)

		# Create flattened data
		self.__create_flattened_data()



	def __parse_output(
		self,
		output_lines
	):
		retdict = {}
		for l in output_lines:
			m=re.match("^\[DATA\] ([^ :]*): (.*)$", l)
			if m != None:
				tag = m.group(1)
				data = m.group(2)
				retdict[tag] = data

		return retdict



	def __load_job_raw_data(
			self,
			jobs_dir = None
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

		if jobs_dir != None:
			jobdirs = glob.glob(jobs_dir+'/job_*')
		else:
			jobdirs = glob.glob('job_*')
			
		# 
		self.__jobs_raw_data = {}
		for jobdir in jobdirs:
			if self.verbosity > 5:
				self.info("Processing '"+jobdir+"'")
			job_data = {}

			if self.verbosity > 5:
				self.info("Loading output 'outout.out'")
			outfile = jobdir+'/output.out'
			with open(outfile) as f:
				content = f.readlines()
				job_data['output'] = self.__parse_output(content)

			if self.verbosity > 5:
				self.info("Loading job generation data 'jobgeneration.pickle'")
			picklefile = jobdir+'/jobgeneration.pickle'
			j = SWEETJobGeneration(dummy_init=True)
			job_data['jobgeneration'] = j.load_attributes_dict(picklefile)

			self.__jobs_raw_data[jobdir] = job_data



	def get_jobs_raw_data(self):
		"""
		Return the raw job information data
		Warning: This storage format is likely to change!
		"""
		return self.__jobs_raw_data



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

		for job_key, job_raw_data in self.__jobs_raw_data.items():
			flattened_job_data = {}

			for key, data in job_raw_data['output'].items():
				key = key.replace('SimulationBenchmarkTimings', 'benchmark_timings')
				key = 'output.'+key.lower()
				flattened_job_data[key] = data

			attr_dict = job_raw_data['jobgeneration']
			for group_name, group_value in attr_dict.items():
				if isinstance(group_value, dict):
					for key, attr in group_value.items():
						key = group_name+"."+key
						flattened_job_data[key] = attr
				elif isinstance(group_value, (int, float, str)):
					flattened_job_data[group_name] = group_value
				else:
					raise Exception("Unknown type in pickled data")

			self.__flattened_jobs_data[job_key] = flattened_job_data

	


	def get_flattened_data(self):
		return self.__flattened_jobs_data


	def create_groups(
		self,
		group_identifiers : list
	):
		"""
		Group together particular job data, e.g. jobs using the same time stepping method

		Parameters:
		-----------
			group_attributes: list
				List with attributes which are identical within jobs of the same group
		"""

		groups = {}

		for jobdir, job_data in self.__flattened_jobs_data.items():
			group_attributes = []
			group_attributes_short = []

			# collect values of group attributes
			for i in group_identifiers:
				group_attributes_short.append(str(job_data[i]))
				group_attributes.append(i+'_'+str(job_data[i]))

			# create string identifier
			group_attributes_id = "_".join(group_attributes)
			group_attributes_short_id = "_".join(group_attributes_short)

			# create new group if it doesn't exist
			if not group_attributes_short_id in groups:
				groups[group_attributes_short_id] = {}

			# append job data to group
			groups[group_attributes_short_id][jobdir] = job_data


		return groups


	def print_data_table(self, data_table):
		for row in data_table:
			print("\t".join([str(d) for d in row]))


	def create_data_plotting(
		self,
		group_identifiers : list,
		primary_key_attribute_name : str,
		data_attribute_name : str,
		placeholder = None
	):
		"""
		Create a data structure which is suitable for plotting

		Parameters:
		-----------
			group_identifiers: list
				List with attributes of jobs which are grouped together
				This also assembles the column name

			primary_attribute_key: str
				Row Key: attribute which should serve as primary key
				(e.g. number of cores, time step size)

			data_attribute_name : str:
				Column data: attribute which should serve as data field

		Return:
		-------
			dict { job_dir : { x_values: list, y_values: list }}
		"""

		
		plot_data = {}

		groups = self.create_groups(group_identifiers)
		for group_id, group_jobs in groups.items():
			x_values = []
			y_values = []
			for jobdir, jobdata in group_jobs.items():
				x = jobdata[primary_key_attribute_name]

				if not data_attribute_name in jobdata:
					y = placeholder
				else:
					y = jobdata[data_attribute_name]

				x_values.append(x)
				y_values.append(y)

			plot_data[group_id] = {
				'x_values': x_values,
				'y_values': y_values
			}

		return plot_data


	def create_data_plotting_float(
		self,
		group_identifiers : list,
		primary_key_attribute_name : str,
		data_attribute_name : str,
		placeholder = None,
		sort_data = True
	):
		data_plotting = self.create_data_plotting(group_identifiers, primary_key_attribute_name, data_attribute_name, placeholder)

		for key, values in data_plotting.items():
			x = [float(i) for i in values['x_values']]
			y = [float(i) for i in values['y_values']]

			if sort_data:
				y = [y for _,y in sorted(zip(x,y))]
				x.sort()

			values['x_values'] = x
			values['y_values'] = y

		return data_plotting
		


	def create_data_table(
		self,
		group_identifiers : list,
		primary_key_attribute_name : str,
		data_attribute_name : str,
		placeholder = None,
		sort_data = True
	):
		"""
		Create a table-like data structure

		Parameters:
		-----------
			group_identifiers: list
				List with attributes of jobs which are grouped together
				This also assembles the column name

			primary_attribute_key: str
				Row Key: attribute which should serve as primary key
				(e.g. number of cores, time step size)

			data_attribute_name : str:
				Column data: attribute which should serve as data field
		"""

		
		groups = self.create_groups(group_identifiers)

		#
		# Determine full set of primary keys in case that primary key is missing somewhere
		#
		row_keys = []
		for group_id, group_jobs in groups.items():
			for jobdir, jobdata in group_jobs.items():
				primary_key = jobdata[primary_key_attribute_name]
				if not primary_key in row_keys:
					row_keys.append(primary_key)
		# Sort the primary keys
		# Not sure if this is always a good idea, but it makes sense for plots with numerical values
		if (sort_data):
			row_keys.sort()

		# get column names
		col_keys = list(groups.keys())

		# get dimensions of table
		ncols = len(groups)
		nrows = len(row_keys)

		# Create table data
		data = [[None for i in range(ncols+1)] for j in range(nrows+1)]

		for group_id, group_jobs in groups.items():

			col_key = group_id
			col_id = col_keys.index(col_key)

			for jobdir, jobdata in group_jobs.items():
				row_key = jobdata[primary_key_attribute_name]
				row_id = row_keys.index(row_key)

				if data[row_id+1][col_id+1] != None:
					self.print_data_table(data)
					raise Exception("Duplicate entry detected! Stopping here")

				data[row_id+1][col_id+1] = jobdata[data_attribute_name]

		data[0][0] = '-'
		#data[0][0] = primary_key_attribute_name+'\\'+data_attribute_name

		# Setup row labels (primary keys)
		for i in range(len(row_keys)):
			data[i+1][0] = row_keys[i]

		# Setup col labels
		for j in range(len(col_keys)):
			data[0][j+1] = col_keys[j]

		if placeholder != None:
			# Replace None fields with placeholder
			for j in range(1,len(row_keys)):
				for i in range(1,len(col_keys)):
					if data[j][i] == None:
						data[j][i] = placeholder

		return data


	def create_data_table_float(
		self,
		group_identifiers : list,
		primary_key_attribute_name : str,
		data_attribute_name : str,
		placeholder = None
	):
		data = self.create_data_table(group_identifiers, primary_key_attribute_name, data_attribute_name)

		# Replace None fields with placeholder
		for j in range(1,len(data)):
			for i in range(1,len(data[0])):
				if data[j][i] == None:
					data[j][i] = placeholder
				else:
					data[j][i] = float(data[j][i])

		return data

if __name__ == "__main__":
	if False:
		j = SWEETPostprocessing(verbosity=0)

		jobs_raw_data = j.get_jobs_raw_data()
		for key, values in jobs_raw_data.items():
			print(key)
			print(values)
			print("")

	#if True:
	if False:
		j = SWEETPostprocessing(verbosity=0)
		d = j.get_flattened_data()

		for jobdir, job_data in d.items():
			print("")
			print("Data for '"+jobdir+"'")
			for key, value in job_data.items():
				print(key+" => "+str(value))

	#if False:
	if True:
		j = SWEETPostprocessing(verbosity=0)
		d = j.get_flattened_data()

		for jobdir, job_data in d.items():
			print("")
			print("Data for '"+jobdir+"'")
			for key, value in job_data.items():
				print(key+" => "+str(value))
			break

	#if True:
	if False:
		j = SWEETPostprocessing(verbosity=0)
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

		for group_id, group_jobs in g.items():
			print(group_id)
			for jobdir, jobdata in group_jobs.items():
				print(" + "+jobdir+"\t"+jobdata['output.benchmark_timings.main_simulationloop'])


	if True:
	#if False:
		j = SWEETPostprocessing(verbosity=0)

		data_table = j.create_data_table(
				['runtime.timestepping_method'],
				'parallelization.num_threads_per_rank',
				'output.benchmark_timings.main_simulationloop'
			)

		print("Data table:")
		j.print_data_table(data_table)
