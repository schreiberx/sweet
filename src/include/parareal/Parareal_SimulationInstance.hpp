/*
 * PararealSimulation.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_HPP_


#include <parareal/Parareal_Data.hpp>
#include "../sweet/plane/PlaneData_Spectral.hpp"


/**
 * Interface descriptions which are required
 * to run Parareal Simulations in SWEET.
 *
 * Note that these interface descriptions have
 * to be implemented by the simulation!
 *
 * These interfaces were ported from the Python implementation.
 */
class Parareal_SimulationInstance
{
public:
	/**
	 * Set the start and end of the coarse time step
	 */
	virtual void sim_set_timeframe(
			double i_timeframe_start,	///< start timestamp of coarse time step
			double i_timeframe_end		///< end time stamp of coarse time step
	) = 0;

	/**
	 * Set the initial data at i_timeframe_start
	 */
	virtual void sim_setup_initial_data(
	) = 0;

	/**
	 * Set simulation data to data given in i_sim_data.
	 * This can be data which is computed by another simulation.
	 * Y^S := i_sim_data
	 */
	virtual void sim_set_data(
			Parareal_Data &i_pararealData
	) = 0;

	/**
	 * Set the MPI communicator to use for simulation purpose
	 * (TODO: not yet implemented since our parallelization-in-space
	 * is done only via OpenMP)
	 */
	virtual void sim_set_mpi_comm(
			int i_mpi_comm
	) = 0;

	/**
	 * compute solution on time slice with fine timestep:
	 * Y^F := F(Y^S)
	 */
	virtual
	void run_timestep_fine() = 0;


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	virtual
	Parareal_Data& get_reference_to_data_timestep_fine() = 0;


	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	virtual
	void run_timestep_coarse() = 0;


	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	virtual
	Parareal_Data& get_reference_to_data_timestep_coarse() = 0;


	/**
	 * Compute the error between the fine and coarse timestepping:
	 * Y^E := Y^F - Y^C
	 */
	virtual
	void compute_difference() = 0;


	/**
	 * Compute the data to be forwarded to the next time step
	 * Y^O := Y^C + Y^E
	 *
	 * Return: If true, the error indicator based on the computed error norm between the
	 * old values and new values
	 */
	virtual
	double compute_output_data(
			bool i_return_convergence
	) = 0;


	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	virtual
	Parareal_Data& get_reference_to_output_data() = 0;

	virtual
	void output_data_console(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	) = 0;

	virtual
	void output_data_file(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	) = 0;

	virtual ~Parareal_SimulationInstance()
	{
	}
};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_HPP_ */
