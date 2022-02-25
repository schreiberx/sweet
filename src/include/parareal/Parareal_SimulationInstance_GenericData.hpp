/*
 * PararealSimulation.hpp
 *
 *  Created on: 25 Feb 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_


#include <parareal/Parareal_Data.hpp>
#include "../sweet/plane/PlaneData.hpp"
#include "../sweet/sphere/SphereData_Spectral.hpp"


/**
 * Interface descriptions which are required
 * to run Parareal Simulations in SWEET.
 *
 * Note that these interface descriptions have
 * to be implemented by the simulation!
 *
 * These interfaces were ported from the Python implementation.
 */
class Parareal_SimulationInstance_GenericData
{
public:

	// Simulation variables
	SimulationVariables simVars;

	// Operators
	PlaneOperators op_plane;
	SphereOperators_SphereData op_sphere;
	SphereOperators_SphereData op_sphere_nodealiasing;

	// Geometry (scalar, plane, sphere)
	std::string geometry = "";

	// Model (ODE, Burgers, SWE)
	std::string model = "";

	// Time slice
	double timeframe_start;
	double timeframe_end;

	// Data containers
	Parareal_Data_GenericData* parareal_data_start = nullptr;
	Parareal_Data_GenericData* parareal_data_fine = nullptr;
	Parareal_Data_GenericData* parareal_data_coarse = nullptr;
	Parareal_Data_GenericData* parareal_data_output = nullptr;
	Parareal_Data_GenericData* parareal_data_error = nullptr;
	Parareal_Data_GenericData* parareal_data_coarse_previous_timestep = nullptr;
	Parareal_Data_GenericData* parareal_data_coarse_previous_time_slice = nullptr;
	Parareal_Data_GenericData* parareal_data_fine_previous_timestep = nullptr;
	Parareal_Data_GenericData* parareal_data_fine_previous_time_slice = nullptr;
#if SWEET_DEBUG
	Parareal_Data_GenericData* parareal_data_fine_exact = nullptr;
#endif

	// Fine and coarse timesteppers
	Generic_TimeSteppers* timeSteppers = nullptr;
	Generic_TimeSteppers* timeSteppersCoarse = nullptr;

	// list of SL schemes
	std::vector<std::string> SL_tsm = {};

	
public:

	// Scalar
	Parareal_SimulationInstance_GenericData(simulationVariables &i_simVars, std::string i_geometry, std::string i_model):
		simVars(i_simVars), geometry(i_geometry), model(i_model)
	{
	};

	// Plane
	Parareal_SimulationInstance_GenericData(simulationVariables &i_simVars, PlaneOperators &i_op_plane, std::string i_geometry, std::string i_model):
		simVars(i_simVars), op_plane(i_op_plane), geometry(i_geometry), model(i_model)
	{
	};

	// Sphere
	Parareal_SimulationInstance_GenericData(simulationVariables &i_simVars, SphereOperators &i_op_sphere, SphereOperators &i_op_sphere_nodealiasing, std::string i_geometry, std::string i_model):
		simVars(i_simVars), op_sphere(i_op_sphere), op_sphere_nodealising(i_op_sphere_nodealising), geometry(i_geometry), model(i_model)
	{
	};



	void setup(simulationVariables &i_simVars, std::string i_geometry, double i_timeframe_start, double i_timeframe_end)
	{

		// Replacing the function sim_set_timeframe
		if (simVars.parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;
		this->timeframe_start = i_timeframe_start;
		this->timeframe_end = i_timeframe_end;


		if (this->geometry == "scalar")
		{
			this->parareal_data_start = new Parareal_GenericData_Scalar;
			this->parareal_data_fine = new Parareal_GenericData_Scalar;
			this->parareal_data_coarse = new Parareal_GenericData_Scalar;
			this->parareal_data_output = new Parareal_GenericData_Scalar;
			this->parareal_data_error = new Parareal_GenericData_Scalar;
			this->parareal_data_coarse_previous_timestep = new Parareal_GenericData_Scalar;
			this->parareal_data_coarse_previous_time_slice = new Parareal_GenericData_Scalar;
			this->parareal_data_fine_previous_timestep = new Parareal_GenericData_Scalar;
			this->parareal_data_fine_previous_time_slice = new Parareal_GenericData_Scalar;
#if SWEET_DEBUG
			this->parareal_data_fine_exact = new Parareal_GenericData_Scalar;
#endif

			this->timeSteppers = new SWE_Plane_TimeSteppers;
			this->timeSteppersCoarse = new SWE_Plane_TimeSteppers;
		}
		else if (this->geometry == "plane")
		{
			this->parareal_data_start = new Parareal_GenericData_PlaneData;
			this->parareal_data_fine = new Parareal_GenericData_PlaneData;
			this->parareal_data_coarse = new Parareal_GenericData_PlaneData;
			this->parareal_data_output = new Parareal_GenericData_PlaneData;
			this->parareal_data_error = new Parareal_GenericData_PlaneData;
			this->parareal_data_coarse_previous_timestep = new Parareal_GenericData_PlaneData;
			this->parareal_data_coarse_previous_time_slice = new Parareal_GenericData_PlaneData;
			this->parareal_data_fine_previous_timestep = new Parareal_GenericData_PlaneData;
			this->parareal_data_fine_previous_time_slice = new Parareal_GenericData_PlaneData;
#if SWEET_DEBUG
			this->parareal_data_fine_exact = new Parareal_GenericData_PlaneData;
#endif

			if (this->model == "Burgers")
			{
				this->timeSteppers = new Burgers_Plane_TimeSteppers;
				this->timeSteppersCoarse = new Burgers_Plane_TimeSteppers;
			}
			else if (this->model == "SWE")
			{
				this->timeSteppers = new SWE_Plane_TimeSteppers;
				this->timeSteppersCoarse = new SWE_Plane_TimeSteppers;
			}
			else
				SWEETError("Unknown model.");

			timeSteppers.setup(
					simVars.disc.timestepping_method,
					simVars.disc.timestepping_order,
					simVars.disc.timestepping_order2,
					op,
					simVars
				);

			timeSteppersCoarse.setup(
					simVars.parareal.coarse_timestepping_method,
					simVars.parareal.coarse_timestepping_order,
					simVars.parareal.coarse_timestepping_order2,
					op,
					simVars
				);

		}
		else if (this->geometry == "sphere")
		{
			this->parareal_data_start = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_fine = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_coarse = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_output = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_error = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_coarse_previous_timestep = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_coarse_previous_time_slice = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_fine_previous_timestep = new Parareal_GenericData_SphereData_Spectral;
			this->parareal_data_fine_previous_time_slice = new Parareal_GenericData_SphereData_Spectral;
#if SWEET_DEBUG
			this->parareal_data_fine_exact = new Parareal_GenericData_SphereData_Spectral;
#endif

			if (this->model == "SWE")
			{
				this->timeSteppers = new SWE_Sphere_TimeSteppers;
				this->timeSteppersCoarse = new SWE_Sphere_TimeSteppers;
			}
			else
				SWEETError("Unknown model.");

			timeSteppers.setup(
					simVars.disc.timestepping_method,
					op,
					simVars
				);

			timeSteppersCoarse.setup(
					simVars.parareal.coarse_timestepping_method,
					op,
					simVars
				);

		}
		else
			SWEETError("Unknown geometry.");

	};



public:
	/**
	 * Check if the time slice contains an integer number of coarse and fine time steŝ
	 */
	void sim_check_timesteps(
			double time_slice_size
	)
	{
		// check if each time slice contains an integer number of fine and coarse time steps
		double eps = 1e-12;
		double mod_coarse = fmod(time_slice_size, simVars.parareal.coarse_timestep_size);
		double mod_fine = fmod(time_slice_size, simVars.timecontrol.current_timestep_size);
                if ( std::abs(mod_coarse) > eps && std::abs(mod_coarse - time_slice_size) > eps  )
			SWEETError("Time slice length must be an integer multiple of the coarse time step!");
                if ( std::abs(mod_fine) > eps && std::abs(mod_fine - time_slice_size) > eps  )
			SWEETError("Time slice length must be an integer multiple of the fine time step!");
	};


	/**
	 * Set the start and end of the coarse time step
	 */
	void sim_set_timeframe(
			double i_timeframe_start,	///< start timestamp of coarse time step
			double i_timeframe_end		///< end time stamp of coarse time step
	){
	};

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
	void sim_set_data(
			Parareal_Data &i_pararealData
	){
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data()" << std::endl;

		// copy to buffers
		parareal_data_start = i_pararealData;
	};

	/**
	 * Set solution of penult coarse timestep of previous time slice
	 */
	void sim_set_data_coarse_previous_time_slice(
			Parareal_Data &i_pararealData
	){
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data_coarse_previous_time_slice()" << std::endl;

		// copy to buffers
		parareal_data_coarse_previous_time_slice = i_pararealData;
	};

	/**
	 * Set solution of penult fine timestep of previous time slice
	 */
	virtual void sim_set_data_fine_previous_time_slice(
			Parareal_Data &i_pararealData
	){
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data_fine_previous_time_slice()" << std::endl;

		// copy to buffers
		parareal_data_fine_previous_time_slice = i_pararealData;
	};

#if SWEET_DEBUG
	/**
	* Store exact solution (full fine simulation) at the end of the time slice
	*/
	virtual void sim_set_data_fine_exact(
			Parareal_Data &i_pararealData
	) = 0;

	/**
	* Check if solution at time k (end of time slice k-1) is exact (= fine) at iteration k
	*/
	virtual void compare_to_fine_exact(
	) = 0;

#endif

	void set_previous_solution()
	{
		if (this->geometry == "scalar")
		else if (this->geometry == "plane")
		else if (this->geometry == "sphere")
	};

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
	void run_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_fine()" << std::endl;

		prog_h_pert = *parareal_data_start.data_arrays[0];
		prog_u = *parareal_data_start.data_arrays[1];
		prog_v = *parareal_data_start.data_arrays[2];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if ( if std::find(this->SL_tsm, simVars.disc.timestepping_method) != this->SL_tsm.end() )
		{
			///PlaneData h_prev = *parareal_data_fine_previous_time_slice.data_arrays[0];
			///PlaneData u_prev = *parareal_data_fine_previous_time_slice.data_arrays[1];
			///PlaneData v_prev = *parareal_data_fine_previous_time_slice.data_arrays[2];

			this->set_previous_solution();
			
			///timeSteppers.master->set_previous_solution(h_prev, u_prev, v_prev);
		}

		while (simVars.timecontrol.current_simulation_time != timeframe_end)
		{
			run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// copy to buffers
		*parareal_data_fine.data_arrays[0] = prog_h_pert;
		*parareal_data_fine.data_arrays[1] = prog_u;
		*parareal_data_fine.data_arrays[2] = prog_v;
	};


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
	 * return the penult solution after the coarse propagation:
	 */
	virtual
	Parareal_Data& get_reference_to_data_timestep_coarse_previous_timestep() = 0;

	/**
	 * return the penult solution after the fine propagation:
	 */
	virtual
	Parareal_Data& get_reference_to_data_timestep_fine_previous_timestep() = 0;


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

	virtual void check_for_nan_parareal() = 0;

	virtual ~Parareal_SimulationInstance()
	{
	}
};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_ */
