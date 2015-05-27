#
# Copyright (c) 2015 University of Exeter
# This file is part of the PyPinT project. For conditions of distribution and
# use, please see the copyright notice in the file 'copyright.txt' at the root
# directory of this package and the copyright notice at
# https://bitbucket.org/schreiberx/pypint
#
#      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
#

from libpypint.PintSimulation_Base import PintSimulation_Base
from libpypint.PintPlot_Base import PintPlot_Base
import matplotlib.pyplot as plt

import math
import copy
import numpy as np
from scipy import interpolate
import sys


class PintSimulation(PintSimulation_Base, PintPlot_Base):
	"""
	Compute solution of 1D SWE PDE
	"""


	"""
	Plot limits for y axis (automatically determined and fixed after first plot)
	"""
	plot_ylim = [None, None]



	"""
	Setup canvas and output filename
    """
	def plot_preprocess(
					self,
					i_filename,
					i_title,
					i_time_interval
	):
		self.plot_filename = i_filename

		plt.clf()
		plt.xlim(0, self.domain_size)

		if self.plot_ylim[0] != None:
			plt.ylim(self.plot_ylim[0], self.plot_ylim[1])

		plt.xlabel('time (x)')
		plt.ylabel('solution (y)')
		plt.title(i_title);#'Iteration '+str(k))



	def plot_start_solution(
				self,
				i_color="red"
	):
		inv_N = 1.0/float(self.N)
		plt.plot([(i+0.5)*inv_N*self.domain_size for i in range(self.N)], self.sim_data_start[0], color=i_color)



	def plot_fine_solution(
				self,
				i_color="red"
	):
		inv_N = 1.0/float(self.N)
		plt.plot([(i+0.5)*inv_N*self.domain_size for i in range(self.N)], self.sim_data_end_fine[0], color=i_color)



	def plot_coarse_solution(
			self,
			i_color="black"
	):
		inv_N = 1.0/float(self.N)
		plt.plot([(i+0.5)*inv_N*self.domain_size for i in range(self.N)], self.sim_data_end_coarse[0], color=i_color)



	"""
	Close filehandle for plotting and set canvas
    """
	def plot_postprocess(
			self
	):
		plt.savefig(self.plot_filename)

		if self.plot_ylim[0] == None:
			(self.plot_ylim[0], self.plot_ylim[1]) = plt.ylim()
		plt.clf()



	"""
	number of grid cell points
	"""
	N = 128

	"""
	domain size in meter
	"""
	domain_size = 500

	"""
	domain size in meter
	"""
	radial_breaking_dam_radius = 50

	"""
	gravitation
	"""
	G = 9.81		

	"""
	Initial height
	"""	
	h0 = 10
	

	"""
	Displacement height
	"""	
	hd = 0.1


	"""
	coarse interpolation method
	"""
#	coarse_interpolation_method = "simple"
	coarse_interpolation_method = "spline"


	"""
	Iteration counter for internal purpose
	"""
#private
	fine_iteration_counter = 0

	"""
	The simulation data:
	For an ODE, it only consists out of a single scalar value
	"""
	sim_data_start = None

	"""
	The simulation after fine timestepping, based on sim_data_start
	"""	
	sim_data_end_fine = None

	"""
	The simulation after coarse timestepping, based on sim_data_start
	"""	
	sim_data_end_coarse = None
	
	"""
	Storage for computing the error between coarse and fine solution
	"""
	sim_data_coarse_error_correction = None

	"""
	The simulation data to be forwarded to next time interval:
	"""	
	sim_data_end_output = None


	"""
	timeframe of coarse time step:
	[start time stamp, end time stamp]
	"""
	timeframe = [-1, -1]


	"""
	number of coarse cells to restrict
	"""
	restrict_num_cells = 2


	"""
	initial values
	"""
	initial_scenario = "square"

	"""
	last coarse timestep size
	"""
	last_coarse_dt = None

	"""
	non-linear or linearized
	"""
	swe_eq_type = "non_linear"


	"""
	time step size for fine time steps
	"""
#private
#	timestep_fine_size = -1


#public
	def __init__(
				self,
				i_argv,
				i_pintOutput = None
	):
		PintPlot_Base.__init__(self, i_argv)

		self.setPintOutput(i_pintOutput)

		for i in range(1, len(i_argv)):
			if i_argv[i]=="res":
				self.N = int(i_argv[i+1]);
				break
			
		for i in range(1, len(i_argv)):
			if i_argv[i]=="rnc":
				self.restrict_num_cells = int(i_argv[i+1]);
				break

		for i in range(1, len(i_argv)):
			if i_argv[i]=="swe_eq_type":
				self.swe_eq_type = i_argv[i+1];
				break

		for i in range(1, len(i_argv)):
			if i_argv[i]=="isc":
				self.initial_scenario = i_argv[i+1];
				break
			
		for i in range(1, len(i_argv)):
			if i_argv[i]=="cim":
				self.coarse_interpolation_method = i_argv[i+1];
				break


		if self.N % self.restrict_num_cells != 0:
			print "N has to be dividable by number of cells involved in restriction with remainder 0"
			sys.exit(-1) 



	"""
	Set the coarse time frame
	"""
	def set_simulation_timeframe(
				self,
				i_timeframe,
				i_coarse_time_frames = 1
	):
		self.pintOutput.outputHeadline("set_simulation_timeframe")
		self.timeframe = i_timeframe[:]
		return



	"""
	Setup the initial values
	"""
	def setup_initial_values(self):

		self.pintOutput.outputHeadline("setup_initial_values")

		self.sim_data_start = np.array(
				[
					[self.h0 for i in range(self.N)],
					[0 for i in range(self.N)]
				]
				, np.double
			)

		cell_size = self.domain_size/float(self.N)

		if self.initial_scenario == "square":
			for i in range(self.N):
				cell_center = cell_size*(i+0.5)
	
				# compute distance to center of domain
				if abs(self.domain_size*0.5 - cell_center) < self.radial_breaking_dam_radius:
					self.sim_data_start[0][i] = self.sim_data_start[0][i]+self.hd

		elif self.initial_scenario == "square2":
			for i in range(self.N):
				cell_center = cell_size*(i+0.5)
	
				# compute distance to center of domain
				if abs(self.domain_size*0.75 - cell_center) < self.radial_breaking_dam_radius:
					self.sim_data_start[0][i] = self.sim_data_start[0][i]+self.hd

		elif self.initial_scenario == "gauss":
			for i in range(self.N):
				cell_center = cell_size*(i+0.5)

				sigma = 0.1	
				mu = 0.5
				x = cell_center/self.domain_size

				dist = 1.0/(sigma*math.sqrt(2.0*math.pi))*math.exp(-(x-mu)**2/(2*sigma*sigma))
				
				dist *= 0.3*self.hd

				self.sim_data_start[0][i] = self.sim_data_start[0][i]+dist

		else:
				assert False, "Unknown scenario"
			

	"""
	Set the simulation data
	"""
	def set_simulation_data(
						self,
						i_sim_data_start
	):
		self.pintOutput.outputHeadline("set_simulation_data")

		self.sim_data_start = copy.deepcopy(i_sim_data_start);
		self.pintOutput.output("  value: "+str(self.sim_data_start))
		return


	def set_mpi_communicator(
					self,
					i_mpi_communicator
	):
		"""dummy"""


	"""
	compute flux for SWE:
	
	U := (h,hu)

	F(U) := (h u, h u^2 + 0.5 g h^2)^T
	
	d/dt U + d/dx F(U) = 0
	"""
#private
	def run_timestep(
				self,
				i_sim_data,
				i_max_dt,
				i_cell_size,
				o_sim_data
	):
		if self.swe_eq_type == "non_linear":
			assert i_sim_data[0][0] >= 0, "Height is equal or less to zero"

			h = i_sim_data[0]
			hu = i_sim_data[1]
			
			q0 = h
			q1 = hu

			u = i_sim_data[1] / i_sim_data[0]

			"compute fluxes in all cells"
			"F(U) := (h u, h u^2 + 0.5 g h^2)^T"
			flux_q0 = hu
			flux_q1 = hu*u + h*h*0.5*self.G
	
			"maximum values of eigenvalues"
			max_eigenvalues = np.absolute(u)+np.sqrt(h*self.G)

		elif self.swe_eq_type == "linearized":
		
			h = i_sim_data[0]
			u = i_sim_data[1]
			
			q0 = h
			q1 = u

			flux_q0 = u*self.h0
			flux_q1 = h*self.G

			max_eigenvalues = np.empty(len(i_sim_data[0]))
			max_eigenvalues.fill(max(self.G, self.h0))

		else:
			assert False, "invalid SWE equation type"
		
		max_speed = i_cell_size/np.max(max_eigenvalues)

		# run with CFL 0.5
		CFL=0.5
		dt_limitation = max_speed * CFL

		dt = min(dt_limitation, i_max_dt)
		
		# important: use this copy operation, otherwise the reference to o_sim_data is lost!
		o_sim_data[:] = i_sim_data[:]

		A = dt/i_cell_size

		local_N = len(i_sim_data[0])
		for i in range(local_N):
			li = local_N-1 if i == 0 else i-1
			mi = i
			ri = 0 if i == local_N-1 else i+1

			if False:
				o_sim_data[0][i] -= A*0.5*(
								-(flux_q0[li]+flux_q0[mi]) - max(max_eigenvalues[li], max_eigenvalues[mi])*(q0[li]-q0[mi])
								+(flux_q0[mi]+flux_q0[ri]) + max(max_eigenvalues[mi], max_eigenvalues[ri])*(q0[mi]-q0[ri])
						)
				o_sim_data[1][i] -= A*0.5*(
								-(flux_q1[li]+flux_q1[mi]) - max(max_eigenvalues[li], max_eigenvalues[mi])*(q1[li]-q1[mi])
								+(flux_q1[mi]+flux_q1[ri]) + max(max_eigenvalues[mi], max_eigenvalues[ri])*(q1[mi]-q1[ri])
						)
			else:
				# LEFT EDGE
				o_sim_data[0][i] -= A*(
								-0.5*(flux_q0[li]+flux_q0[mi]) - 0.5*max(max_eigenvalues[li], max_eigenvalues[mi])*(q0[li]-q0[mi])
						)
				o_sim_data[1][i] -= A*(
								-0.5*(flux_q1[li]+flux_q1[mi]) - 0.5*max(max_eigenvalues[li], max_eigenvalues[mi])*(q1[li]-q1[mi])
						)
	
				# RIGHT EDGE: TODO: unify with left edge
				o_sim_data[0][i] -= A*(
								0.5*(flux_q0[mi]+flux_q0[ri]) + 0.5*max(max_eigenvalues[mi], max_eigenvalues[ri])*(q0[mi]-q0[ri])
						)
				o_sim_data[1][i] -= A*(
								0.5*(flux_q1[mi]+flux_q1[ri]) + 0.5*max(max_eigenvalues[mi], max_eigenvalues[ri])*(q1[mi]-q1[ri])
						)
		return dt



	"""
	For sake of simplicity, we use an explicit Euler timestepping
	method here.
	This can be replaced by any others since it's transparent to this interface
	"""
	def run_timestep_fine(self):
		self.pintOutput.outputHeadline("run_timestep_fine with iteration nr. "+str(self.fine_iteration_counter))
		self.fine_iteration_counter = self.fine_iteration_counter+1

#		self.sim_data_end_fine = np.zeros((2, self.N), np.double)
		self.sim_data_end_fine = np.copy(self.sim_data_start)

		assert self.timeframe[0] != -1, "timeframe not set (set_simulation_data not executed)"

		t = self.timeframe[0]
		while t < self.timeframe[1]:
			# copy initial values to temporary start buffer
			tmp_start = np.copy(self.sim_data_end_fine)

			"compute the maximum allowed time step size"
			max_dt = self.timeframe[1]-t

			real_dt = self.run_timestep(
							tmp_start,
							max_dt,
							self.domain_size/float(self.N),
							self.sim_data_end_fine
						)

			t += real_dt

#		print str(self.last_coarse_dt)+" "+str(real_dt)
		if self.last_coarse_dt != None:
			if self.last_coarse_dt*0.5 <= real_dt:
				print "WARNING: Coarse time step size ("+str(self.last_coarse_dt)+") of similar magnitude as fine time step size ("+str(real_dt)+")!"
#				self.last_coarse_dt

		self.pintOutput.output("  start: "+str(self.sim_data_end_fine))


		"""
		if self.output_plots:
			self.debug_fine_timestep_data_t = []
			self.debug_fine_timestep_data_y = []

		t = self.timeframe[0]
		while t < self.timeframe[1]:
			if self.output_plots:
				self.debug_fine_timestep_data_t.append(t)
				self.debug_fine_timestep_data_y.append(self.sim_data_end_fine)

			"limit fine timestep size to avoid overshooting the last timestep"
			dt = min(self.timestep_fine_size, self.timeframe[1]-t)

			"simulation data of next timestep"
			self.sim_data_end_fine += dt*self.eval_df_dt(t, self.sim_data_end_fine)

			t += dt

		if self.output_plots:
			self.debug_fine_timestep_data_t.append(t)
			self.debug_fine_timestep_data_y.append(self.sim_data_end_fine)

		self.pintOutput.output("  end: "+str(self.sim_data_end_fine))
		"""



	def get_data_timestep_fine(self):
		self.pintOutput.outputHeadline("get_data_timestep_fine")
		self.pintOutput.output("  value: "+str(self.sim_data_end_fine))
		return self.sim_data_end_fine



	"""
	Here, we give it a try with the explicit Euler timestepping method.
	Let's hope that we are still converging...
	"""
	def run_timestep_coarse(self):
		self.pintOutput.outputHeadline("run_timestep_coarse")

		"""
		RESTRICT fine solution to allow larger TS
		"""
		num_coarse_cells = self.N/self.restrict_num_cells

		cell_size_fine = self.domain_size/float(self.N)
		cell_size_coarse = self.domain_size/float(num_coarse_cells)

		coarse_solution = np.zeros((2, num_coarse_cells), np.double)

		if self.coarse_interpolation_method == "simple":
			for i in range(num_coarse_cells):
				for j in range(self.restrict_num_cells):
					coarse_solution[0][i] += self.sim_data_start[0][i*self.restrict_num_cells+j]
					coarse_solution[1][i] += self.sim_data_start[1][i*self.restrict_num_cells+j]
					
			coarse_solution /= float(self.restrict_num_cells)

		elif self.coarse_interpolation_method == "linear":
			x = np.arange(cell_size_fine*0.5, self.domain_size, cell_size_fine)
			xcoarse = np.arange(cell_size_coarse*0.5, self.domain_size, cell_size_coarse)

			for i in [0, 1]:
				interf = interpolate.interp1d(x, self.sim_data_start[i])
				coarse_solution[i] = interf(xcoarse)

		elif self.coarse_interpolation_method == "spline":
			x = np.arange(cell_size_fine*0.5, self.domain_size, cell_size_fine)
			xcoarse = np.arange(cell_size_coarse*0.5, self.domain_size, cell_size_coarse)

			for i in [0, 1]:
				tck = interpolate.splrep(x, self.sim_data_start[i], s=0)
				coarse_solution[i] = interpolate.splev(xcoarse, tck, der=0)

		else:
			assert False, "Unknown interpolation method "+str(self.coarse_interpolation_method)


		"""
		COARSE TIME STEPPING on coarse_solution
		"""
		assert self.timeframe[0] != -1, "timeframe not set (set_simulation_data not executed)"

		t = self.timeframe[0]
		while t < self.timeframe[1]:

			"compute the maximum allowed time step size"
			max_dt = self.timeframe[1]-t

			tmp_end = np.zeros((2, num_coarse_cells), np.double)

			"run timestep and return used time step size"
			real_dt = self.run_timestep(
							coarse_solution,
							max_dt,
							self.domain_size/float(num_coarse_cells),
							tmp_end
						)

			coarse_solution = np.copy(tmp_end)
			t += real_dt
			
		self.last_coarse_dt = real_dt

		"""
		INTERPOLATE coarse solution to fine solution
		
		IMPORTANT: Care about periodic boundary conditions!!!
		"""
		self.sim_data_end_coarse = np.zeros((2, self.N), np.double)

		"extend coarse sampling points by periodic boundary conditions"
		xcoarse_ext = [-cell_size_coarse*0.5]
		xcoarse_ext.extend(xcoarse)
		xcoarse_ext.extend([self.domain_size+cell_size_coarse*0.5])

		"extend coarse solution to periodic boundary conditions"
		coarse_solution_ext = np.zeros((2,num_coarse_cells+2), np.double)
		for i in [0, 1]:
			coarse_solution_ext[i][0] = coarse_solution[i][num_coarse_cells-1]
			coarse_solution_ext[i][1:num_coarse_cells+1] = coarse_solution[i]
			coarse_solution_ext[i][num_coarse_cells+1] = coarse_solution[i][0]

		if self.coarse_interpolation_method == "simple":

			for i in range(num_coarse_cells):
				for j in range(self.restrict_num_cells):
					self.sim_data_end_coarse[0][i*self.restrict_num_cells+j] = coarse_solution[0][i]
					self.sim_data_end_coarse[1][i*self.restrict_num_cells+j] = coarse_solution[1][i]

		elif self.coarse_interpolation_method == "linear":
			for i in [0, 1]:
				interf = interpolate.interp1d(xcoarse_ext, coarse_solution_ext[i])
				self.sim_data_end_coarse[i] = interf(x)

		elif self.coarse_interpolation_method == "spline":
			for i in [0, 1]:
				tck = interpolate.splrep(xcoarse_ext, coarse_solution_ext[i], s=0)
				self.sim_data_end_coarse[i] = interpolate.splev(x, tck, der=0)

		else:
			assert False, "Unknown interpolation method "+str(self.coarse_interpolation_method)


	def get_data_timestep_coarse(self):
		self.pintOutput.outputHeadline("get_data_timestep_coarse")
		self.pintOutput.output("  value: "+str(self.sim_data_end_coarse))
		
		return self.sim_data_end_coarse



	"""
	compute the difference between the fine and coarse timestep solution
	"""
	def compute_coarse_timestep_error(
					self,
	):
		self.pintOutput.outputHeadline("compute_coarse_timestep_error")

		self.sim_data_coarse_error_correction = self.sim_data_end_fine - self.sim_data_end_coarse

		self.pintOutput.output("   sim_data_coarse_error_correction: "+str(self.sim_data_coarse_error_correction))



	"""
	update the currently stored solution with the propagated sim_data_coarse_error_correction
	"""
	def compute_output_data(
			self,
			i_compute_convergence_indicator = False
	):
		self.pintOutput.outputHeadline("compute_output_data")
		
		self.pintOutput.outputHeadline("compute_output_data")

		if not i_compute_convergence_indicator:
			# fix for next time step
			self.sim_data_end_output = self.sim_data_end_coarse + self.sim_data_coarse_error_correction
			return None

		if self.sim_data_end_output == None:
			"if there's no output, return infinity to request one more iteration for convergence"
			self.sim_data_end_output = self.sim_data_end_coarse + self.sim_data_coarse_error_correction
			return float('inf')

		new_output = self.sim_data_end_coarse + self.sim_data_coarse_error_correction
		convergence_indicator = self.compute_norm_on_solutions(new_output, self.sim_data_end_output)
		self.sim_data_end_output = new_output
		return convergence_indicator



	"""
	return data to be forwarded as input data for next coarse time interval
	"""
	def get_output_data(self):
		return self.sim_data_end_output



	"""
	return an sim_data_coarse_error_correction estimation
	"""
	def return_error_estimation(self):
		self.pintOutput.outputHeadline("return_error_estimation")

		return None



	def compute_norm_on_solutions(
				self,
				solution_a,
				solution_b
	):
		value = 0.0

		for i in range(self.N):
			value = value + abs(solution_a[0][i] - solution_b[0][i])		

#		for i in range(self.N):
#			e = abs(solution_a[0][i] - solution_b[0][i])
#			value = value + e*e

		return value / float(self.N)

		