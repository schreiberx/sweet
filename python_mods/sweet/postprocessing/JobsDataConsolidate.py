#! /usr/bin/env python3

from SWEET import *
from sweet.postprocessing.JobData import *
from sweet.postprocessing.JobsData import *
import glob
import os
import re
import sys




class JobsDataConsolidate(InfoError):
	def __init__(
			self,
			verbosity = 0
	):
		self.verbosity = verbosity

	def create_groups(
		self,
		jobs_data,
		group_identifiers : list,
	):
		"""
		Group together particular job data, e.g. jobs using the same time stepping method

		Parameters:
		-----------
			group_attributes: list
				List with attributes which are identical within jobs of the same group
		"""

		groups = {}

		for jobdir, job_data in jobs_data.__flattened_jobs_data.items():
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

	def write_data_table(self, data_table, filename):
		with open(filename, 'w') as f:
			for row in data_table:
				f.write("\t".join([str(d) for d in row])+"\n")

	def create_data_plotting(
		self,
		jobs_data,
		group_identifiers : list,
		primary_key_attribute_name : str,
		data_attribute_name : str,
		placeholder = None,
		data_filter = None
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

		groups = self.create_groups(jobs_data, group_identifiers)
		for group_id in sorted(groups):
			group_jobs = groups[group_id]

			x_values = []
			y_values = []
			for jobdir, jobdata in group_jobs.items():
				if not primary_key_attribute_name in jobdata:
					print("")
					print("WARNING: No data for attribute "+primary_key_attribute_name+" found")
					print("WARNING: Job: "+jobdir)
					continue

				x = jobdata[primary_key_attribute_name]

				if not data_attribute_name in jobdata:
					y = placeholder
				else:
					y = jobdata[data_attribute_name]

					if data_filter != None:
						if data_filter(group_id, x, y):
							continue

				x_values.append(x)
				y_values.append(y)

			plot_data[group_id] = {
				'x_values': x_values,
				'y_values': y_values
			}

		return plot_data



	def print_data_plotting(self, data_plotting):
		for group_key in sorted(data_plotting):
			print("Group '"+group_key+"'")
			group_data = data_plotting[group_key]

			for x, y in zip(group_data['x_values'], group_data['y_values']):
				print("	"+str(x)+" -> "+str(y))
			print("")



	def write_data_plotting(self, data_plotting, filename):
		with open(filename, 'w') as f:
			for group_key in sorted(data_plotting):
				f.write("Group '"+group_key+"'\n")
				group_data = data_plotting[group_key]

				f.write("x_values\ty_values\n")
				for x, y in zip(group_data['x_values'], group_data['y_values']):
					f.write(str(x)+"\t"+str(y)+"\n")
				f.write("\n")



	def create_data_plotting_float(
		self,
		jobs_data,
		group_identifiers : list,
		primary_key_attribute_name : str,
		data_attribute_name : str,
		placeholder = None,
		sort_data = True,
		data_filter = None

	):
		data_plotting = self.create_data_plotting(jobs_data, group_identifiers, primary_key_attribute_name, data_attribute_name, placeholder, data_filter = data_filter)

		for key, values in data_plotting.items():
			x = []
			y = []
			for (i, j) in zip(values['x_values'], values['y_values']):
				if j != None:
					x.append(float(i))
					y.append(float(j))

			if sort_data:
				y = [y for _,y in sorted(zip(x,y))]
				x.sort()

			values['x_values'] = x
			values['y_values'] = y

		return data_plotting
		


	def create_data_table(
		self,
		jobs_data,
		group_identifiers : list,
		primary_key_attribute_name : str,
		data_attribute_name : str,
		placeholder = None,
		sort_data = True,
		data_filter = None
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

		
		groups = self.create_groups(jobs_data, group_identifiers)

		#
		# Determine full set of primary keys in case that primary key is missing somewhere
		#
		row_keys = []
		for group_id, group_jobs in groups.items():
			for jobdir, jobdata in group_jobs.items():
				if not primary_key_attribute_name in jobdata:
					print("")
					print("WARNING: No data for attribute "+primary_key_attribute_name+" found")
					print("WARNING: Job: "+jobdir)
					continue

				primary_key = jobdata[primary_key_attribute_name]
				if not primary_key in row_keys:
					row_keys.append(primary_key)
		# Sort the primary keys
		# Not sure if this is always a good idea, but it makes sense for plots with numerical values
		if (sort_data):
			row_keys.sort()

		# get column names
		col_keys = sorted(list(groups.keys()))

		# get dimensions of table
		ncols = len(groups)
		nrows = len(row_keys)

		# Create table data
		data = [[None for i in range(ncols+1)] for j in range(nrows+1)]

		for group_id, group_jobs in groups.items():
			col_key = group_id
			col_id = col_keys.index(col_key)

			for jobdir, jobdata in group_jobs.items():
				if self.verbosity > 5:
					print("Job: "+jobdir)

				if not primary_key_attribute_name in jobdata:
					print("")
					print("WARNING: No data for attribute "+primary_key_attribute_name+" found")
					print("WARNING: Job: "+jobdir)
					continue

				row_key = jobdata[primary_key_attribute_name]
				row_id = row_keys.index(row_key)

				if data[row_id+1][col_id+1] != None:
					self.print_data_table(data)
					print("")
					print("ERROR: Duplicate entry detected")
					print("ERROR: This typically happens if either")
					print("ERROR:  a) Groups are colliding")
					print("ERROR:  b) axis variables are incorrect")
					print("")
					raise Exception("Duplicate entry!")

				if not data_attribute_name in jobdata:
					print("WARNING: attribute "+data_attribute_name+" not found")
					print("Job directory: "+jobdir)
					#for key, value in jobdata.items():
					#	print(" + "+key+": "+str(value))

					# Ignore missing data, will be filled in by placeholder :-)
					#raise Exception("attribute '"+data_attribute_name+"' not found")

				else:
					x = row_key
					y = jobdata[data_attribute_name]
					if data_filter != None:
						if data_filter(group_id, x, y):
							continue

					data[row_id+1][col_id+1] = y

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
		placeholder = None,
		data_filter = None
	):
		data = self.create_data_table(group_identifiers, primary_key_attribute_name, data_attribute_name, data_filter=data_filter)

		# Replace None fields with placeholder
		for j in range(1,len(data)):
			for i in range(1,len(data[0])):
				if data[j][i] == None:
					data[j][i] = placeholder
				else:
					data[j][i] = float(data[j][i])

		return data


