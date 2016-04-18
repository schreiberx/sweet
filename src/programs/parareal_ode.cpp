/*
 * parareal_ode.cpp
 *
 *  Created on: 18 Apr 2016
 *      Author: martin
 */

#if !SWEET_PARAREAL
#	error "Parareal not activated"
#endif

#include <stdlib.h>
#include <cmath>

#include <parareal/Parareal.hpp>
#include <parareal/Parareal_Data.hpp>
#include <parareal/Parareal_Data_Scalar.hpp>
#include <parareal/Parareal_Controller_Serial.hpp>

#include <sweet/SimulationVariables.hpp>

SimulationVariables simVars;

double param_parareal_fine_dt = -1;



/*
 * ODE implementation of the Parareal algorithm to test implementations.
 *
 * Usage of the program:
 * --parareal-fine-dt=0.001 --parareal-enabled=1 --parareal-coarse-slices=10 -t 5 --parareal-convergence-threshold=0.0001
 *
 */


class SimulationInstance	:
		public Parareal_SimulationInstance
{
	Parareal_Data_Scalar parareal_data_start;
	Parareal_Data_Scalar parareal_data_fine;
	Parareal_Data_Scalar parareal_data_coarse;
	Parareal_Data_Scalar parareal_data_output;
	Parareal_Data_Scalar parareal_data_error;


	double timeframe_start = -1;
	double timeframe_end = -1;

	bool output_data_valid = false;



	/**
	 * Set the start and end of the coarse time step
	 */
	void sim_set_timeframe(
			double i_timeframe_start,	///< start timestamp of coarse time step
			double i_timeframe_end		///< end time stamp of coarse time step
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;

		timeframe_start = i_timeframe_start;
		timeframe_end = i_timeframe_end;
	}



	/**
	 * Set the initial data at i_timeframe_start
	 */
	void sim_setup_initial_data(
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_setup_initial_data()" << std::endl;

		parareal_data_start.data = 0.123;
	}



	/**
	 * Set simulation data to data given in i_sim_data.
	 * This can be data which is computed by another simulation.
	 * Y^S := i_sim_data
	 */
	void sim_set_data(
			Parareal_Data &i_pararealData
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data()" << std::endl;

		// copy to buffers
		parareal_data_start = i_pararealData;
	}

	/**
	 * Set the MPI communicator to use for simulation purpose
	 * (TODO: not yet implemented since our parallelization-in-space
	 * is done only via OpenMP)
	 */
	void sim_set_mpi_comm(
			int i_mpi_comm
	)
	{
		// NOTHING TO DO HERE
	}



	/**
	 * Function to be integrated parallel in time
	 */
	double f_dt(double y, double t)
	{
		return std::sin(t) + std::sin(y)*0.1;
	}



	/**
	 * compute solution on time slice with fine timestep:
	 * Y^F := F(Y^S)
	 */
	void run_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_fine()" << std::endl;

		double y = parareal_data_start.data;
		double t = timeframe_start;

		while (t != timeframe_end)
		{
			double dt = std::min(param_parareal_fine_dt, timeframe_end-t);

			if (simVars.misc.verbosity > 5)
				std::cout << t << ": " << y << std::endl;

			y += f_dt(y, t)*dt;
			t += dt;

			assert(t <= timeframe_end);
		}

		if (simVars.misc.verbosity > 5)
			std::cout << t << ": " << y << std::endl;

		parareal_data_fine.data = y;
	}


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_Data& get_data_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_data_timestep_fine()" << std::endl;

		return parareal_data_fine;
	}


	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	void run_timestep_coarse()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_coarse()" << std::endl;

		double y = parareal_data_start.data;
		double dt = timeframe_end - timeframe_start;

		if (simVars.misc.verbosity > 5)
			std::cout << timeframe_start << ": " << y << std::endl;

		y += dt*f_dt(y, timeframe_start);

		if (simVars.misc.verbosity > 5)
			std::cout << timeframe_end << ": " << y << std::endl;

		parareal_data_coarse.data = y;
	}



	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_Data& get_data_timestep_coarse()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_data_timestep_coarse()" << std::endl;

		return parareal_data_coarse;
	}



	/**
	 * Compute the error between the fine and coarse timestepping:
	 * Y^E := Y^F - Y^C
	 */
	void compute_difference()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "compute_difference()" << std::endl;

		parareal_data_error.data = parareal_data_fine.data - parareal_data_coarse.data;
	}



	/**
	 * Compute the data to be forwarded to the next time step
	 * Y^O := Y^C + Y^E
	 *
	 * Return: Error indicator based on the computed error norm between the
	 * old values and new values
	 */
	double compute_output_data(
			bool i_compute_convergence_test
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "compute_output_data()" << std::endl;

		double convergence = -1;

		if (!i_compute_convergence_test || !output_data_valid)
		{
			parareal_data_output.data = parareal_data_coarse.data + parareal_data_error.data;

			output_data_valid = true;
			return convergence;
		}


		for (int k = 0; k < 3; k++)
		{
			double tmp = parareal_data_coarse.data + parareal_data_error.data;

			convergence = std::max(
					convergence,
					std::abs(parareal_data_output.data-tmp)
				);

			parareal_data_output.data = tmp;
		}

		std::cout << "                           " << parareal_data_output.data << std::endl;

		output_data_valid = true;
		return convergence;
	}



	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_Data& get_output_data()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_output_data()" << std::endl;

		return parareal_data_output;
	}


	void output_data_console(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
		std::cout << ((Parareal_Data_Scalar&)i_data).data << std::endl;
	}


	void output_data_file(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
//		std::cout << ((Parareal_Data_Scalar&)i_data).data << std::endl;
	}
};


int main(int i_argc, char *i_argv[])
{
	const char *bogus_var_names[] = {
		"parareal-fine-dt",
		nullptr
	};

	simVars.bogus.var[0] = -1;

	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << "Unknown option detected" << std::endl;
		exit(-1);
	}

	param_parareal_fine_dt = simVars.bogus.var[0];

	if (param_parareal_fine_dt <= 0)
	{
		std::cout << "Specify fine time step size via --parareal-fine-dt=[value]" << std::endl;
		return -1;
	}


	if (!simVars.parareal.enabled)
	{
		std::cout << "Activate parareal mode via --parareal0enable=1" << std::endl;
		return -1;
	}

	/*
	 * Allocate parareal controller and provide class
	 * which implement the parareal features
	 */

	Parareal_Controller_Serial<SimulationInstance> parareal_Controller_Serial;

	// setup controller. This initializes several simulation instances
	parareal_Controller_Serial.setup(&simVars.parareal);


	// execute the simulation
	parareal_Controller_Serial.run();

	return 0;
}
