/*
 * PararealSimulation.hpp
 *
 *  Created on: 25 Feb 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_


#include <parareal/Parareal_GenericData.hpp>
#include <parareal/Parareal_GenericData_Scalar.hpp>
#include <parareal/Parareal_GenericData_PlaneData_Spectral.hpp>
#include <parareal/Parareal_GenericData_SphereData_Spectral.hpp>

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>

#include "../../programs/swe_plane_timeintegrators/SWE_Plane_TimeSteppers.hpp"
#include "../../programs/swe_sphere_timeintegrators/SWE_Sphere_TimeSteppers.hpp"
#include "../../programs/burgers_timeintegrators/Burgers_Plane_TimeSteppers.hpp"


/**
 * Interface descriptions which are required
 * to run Parareal Simulations in SWEET.
 *
 * Note that these interface descriptions have
 * to be implemented by the simulation!
 *
 * These interfaces were ported from the Python implementation.
 */
template <class t_tsmType, template <int N> class t_dataType2, int N>
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
	Parareal_GenericData* parareal_data_start = nullptr;
	Parareal_GenericData* parareal_data_fine = nullptr;
	Parareal_GenericData* parareal_data_coarse = nullptr;
	Parareal_GenericData* parareal_data_output = nullptr;
	Parareal_GenericData* parareal_data_error = nullptr;
	Parareal_GenericData* parareal_data_coarse_previous_timestep = nullptr;
	Parareal_GenericData* parareal_data_coarse_previous_time_slice = nullptr;
	Parareal_GenericData* parareal_data_fine_previous_timestep = nullptr;
	Parareal_GenericData* parareal_data_fine_previous_time_slice = nullptr;
#if SWEET_DEBUG
	Parareal_GenericData* parareal_data_fine_exact = nullptr;
#endif

	// Fine and coarse timesteppers
	t_tsmType* timeSteppersFine = nullptr;
	t_tsmType* timeSteppersCoarse = nullptr;

	// list of SL schemes
	std::vector<std::string> SL_tsm = {};

	
public:

	Parareal_SimulationInstance_GenericData()
	{
	};

///	// Scalar
///	Parareal_SimulationInstance_GenericData(simulationVariables &i_simVars, std::string i_geometry, std::string i_model):
///		simVars(i_simVars), geometry(i_geometry), model(i_model)
///	{
///	};
///
///	// Plane
///	Parareal_SimulationInstance_GenericData(simulationVariables &i_simVars, PlaneOperators &i_op_plane, std::string i_geometry, std::string i_model):
///		simVars(i_simVars), op_plane(i_op_plane), geometry(i_geometry), model(i_model)
///	{
///	};
///
///	// Sphere
///	Parareal_SimulationInstance_GenericData(simulationVariables &i_simVars, SphereOperators &i_op_sphere, SphereOperators &i_op_sphere_nodealiasing, std::string i_geometry, std::string i_model):
///		simVars(i_simVars), op_sphere(i_op_sphere), op_sphere_nodealising(i_op_sphere_nodealising), geometry(i_geometry), model(i_model)
///	{
///	};

	// Plane
	void setup(SimulationVariables &i_simVars, std::string i_geometry, std::string i_model,
		   PlaneOperators &i_op_plane,
		   t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{
		this->setup(i_simVars, i_geometry, i_model, i_timeSteppersFine, i_timeSteppersCoarse);
		this->op_plane = i_op_plane;
	}

	// Sphere
	void setup(SimulationVariables &i_simVars, std::string i_geometry, std::string i_model,
		   SphereOperators_SphereData &i_op_sphere, SphereOperators_SphereData &i_op_sphere_nodealiasing,
		   t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{

		this->setup(i_simVars, i_geometry, i_model, i_timeSteppersFine, i_timeSteppersCoarse);
		this->op_sphere = i_op_sphere;
		this->op_sphere_nodealiasing = i_op_sphere_nodealiasing;
	}


	void setup(SimulationVariables &i_simVars, std::string i_geometry, std::string i_model,
			t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{

		this->timeSteppersFine = i_timeSteppersFine;
		this->timeSteppersCoarse = i_timeSteppersCoarse;

///		int aux_nb_fields;
///		if ( ! i_model.compare("simple"))
///			aux_nb_fields = 1;
///		else if ( ! i_model.compare("burgers"))
///			aux_nb_fields = 2;
///		else if ( ! i_model.compare("swe"))
///			aux_nb_fields = 3;
///		else
///			SWEETError("Unknown model");
///		static const int nb_fields(aux_nb_fields);

			this->parareal_data_start = new t_dataType2<N>;
			this->parareal_data_fine = new t_dataType2<N>;
			this->parareal_data_coarse = new t_dataType2<N>;
			this->parareal_data_output = new t_dataType2<N>;
			this->parareal_data_error = new t_dataType2<N>;
			this->parareal_data_coarse_previous_timestep = new t_dataType2<N>;
			this->parareal_data_coarse_previous_time_slice = new t_dataType2<N>;
			this->parareal_data_fine_previous_timestep = new t_dataType2<N>;
			this->parareal_data_fine_previous_time_slice = new t_dataType2<N>;
#if SWEET_DEBUG
			this->parareal_data_fine_exact = new t_dataType2<N>;
#endif

		///	this->timeSteppersFine = new t_tsmType;
		///	this->timeSteppersCoarse = new t_tsmType;



		///const int nb_fields = 3;
		if (this->geometry == "scalar")
		{
//////////			this->parareal_data_start = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_fine = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_coarse = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_output = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_error = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_coarse_previous_timestep = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_coarse_previous_time_slice = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_fine_previous_timestep = new Parareal_GenericData_Scalar<N>;
//////////			this->parareal_data_fine_previous_time_slice = new Parareal_GenericData_Scalar<N>;
//////////#if SWEET_DEBUG
//////////			this->parareal_data_fine_exact = new Parareal_GenericData_Scalar<N>;
//////////#endif
//////////
//////////			this->timeSteppersFine = new SWE_Plane_TimeSteppers;
//////////			this->timeSteppersCoarse = new SWE_Plane_TimeSteppers;
		}
		else if (this->geometry == "plane")
		{
////////////			this->parareal_data_start = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_fine = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_coarse = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_output = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_error = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_coarse_previous_timestep = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_coarse_previous_time_slice = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_fine_previous_timestep = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////			this->parareal_data_fine_previous_time_slice = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////#if SWEET_DEBUG
////////////			this->parareal_data_fine_exact = new Parareal_GenericData_PlaneData_Spectral<N>;
////////////#endif
////////////
////////////			if (this->model == "Burgers")
////////////			{
////////////				this->timeSteppersFine = new Burgers_Plane_TimeSteppers;
////////////				this->timeSteppersCoarse = new Burgers_Plane_TimeSteppers;
////////////			}
////////////			else if (this->model == "SWE")
////////////			{
////////////				this->timeSteppers = new SWE_Plane_TimeSteppers;
////////////				this->timeSteppersCoarse = new SWE_Plane_TimeSteppers;
////////////			}
////////////			else
				SWEETError("Unknown model.");

///			timeSteppersFine->setup(
///					simVars.disc.timestepping_method,
///					simVars.disc.timestepping_order,
///					simVars.disc.timestepping_order2,
///					op_plane,
///					simVars
///				);
///
///			timeSteppersCoarse->setup(
///					simVars.parareal.coarse_timestepping_method,
///					simVars.parareal.coarse_timestepping_order,
///					simVars.parareal.coarse_timestepping_order2,
///					op_plane,
///					simVars
///				);

		}
		else if (this->geometry == "sphere")
		{
//////////////			this->parareal_data_start = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_fine = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_coarse = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_output = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_error = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_coarse_previous_timestep = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_coarse_previous_time_slice = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_fine_previous_timestep = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////			this->parareal_data_fine_previous_time_slice = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////#if SWEET_DEBUG
//////////////			this->parareal_data_fine_exact = new Parareal_GenericData_SphereData_Spectral<N>;
//////////////#endif
//////////////
//////////////			if (this->model == "SWE")
//////////////			{
//////////////				this->timeSteppersFine = new SWE_Sphere_TimeSteppers;
//////////////				this->timeSteppersCoarse = new SWE_Sphere_TimeSteppers;
//////////////			}
//////////////			else
//////////////				SWEETError("Unknown model.");

///			timeSteppersFine->setup(
///					simVars.disc.timestepping_method,
///					op_sphere,
///					simVars
///				);
///
///			timeSteppersCoarse->setup(
///					simVars.parareal.coarse_timestepping_method,
///					op_sphere,
///					simVars
///				);

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
		// Replacing the function sim_set_timeframe
		if (simVars.parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;
		this->timeframe_start = i_timeframe_start;
		this->timeframe_end = i_timeframe_end;
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
			Parareal_GenericData &i_pararealData
	){
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data()" << std::endl;

		// copy to buffers
		*parareal_data_start = i_pararealData;
	};

	/**
	 * Set solution of penult coarse timestep of previous time slice
	 */
	void sim_set_data_coarse_previous_time_slice(
			Parareal_GenericData &i_pararealData
	){
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data_coarse_previous_time_slice()" << std::endl;

		// copy to buffers
		*parareal_data_coarse_previous_time_slice = i_pararealData;
	};

	/**
	 * Set solution of penult fine timestep of previous time slice
	 */
	virtual void sim_set_data_fine_previous_time_slice(
			Parareal_GenericData &i_pararealData
	){
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data_fine_previous_time_slice()" << std::endl;

		// copy to buffers
		*parareal_data_fine_previous_time_slice = i_pararealData;
	};

#if SWEET_DEBUG
	/**
	* Store exact solution (full fine simulation) at the end of the time slice
	*/
	virtual void sim_set_data_fine_exact(
			Parareal_GenericData &i_pararealData
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
		{
		}
		else if (this->geometry == "plane")
		{
		}
		else if (this->geometry == "sphere")
		{
		}
	};

	/**
	 * Set the MPI communicator to use for simulation purpose
	 * (TODO: not yet implemented since our parallelization-in-space
	 * is done only via OpenMP)
	 */
	virtual void sim_set_mpi_comm(
			int i_mpi_comm
	) = 0;


	void run_timestep(
			Parareal_GenericData &i_data,
			std::string tsm_level
	)
	{
		if ( ! (tsm_level == "fine" || tsm_level == "coarse"))
			SWEETError("Wrong tsm_level (should be 'fine' or 'coarse')");


		if (this->geometry == "scalar")
		{
		}
		else if (this->geometry == "plane")
		{
			if (this->model == "burgers")
			{
			}
			else if (this->model == "swe")
			{

				PlaneData prog_h_pert = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
				PlaneData prog_u = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[1]);
				PlaneData prog_v = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[2]);

				if (tsm_level == "fine")
					timeSteppersFine->master->run_timestep(
								prog_h_pert, prog_u, prog_v,
								simVars.timecontrol.current_timestep_size,
								simVars.timecontrol.current_simulation_time
							);
				else if (tsm_level == "coarse")
					timeSteppersCoarse->master->run_timestep(
								prog_h_pert, prog_u, prog_v,
								simVars.timecontrol.current_timestep_size,
								simVars.timecontrol.current_simulation_time
							);

				// copy to buffers
				*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = prog_h_pert;
				*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = prog_u;
				*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[2]) = prog_v;

			}
			else
				SWEETError("Unknown model for this geometry");
		}
		else if (this->geometry == "sphere")
		{
			if (this->model == "swe")
			{
			}
			else
				SWEETError("Unknown model for this geometry");
		}
		else
			SWEETError("Unknown geometry");

	}



	/**
	 * compute solution on time slice with fine timestep:
	 * Y^F := F(Y^S)
	 */
	void run_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_fine()" << std::endl;

		//prog_h_pert = *parareal_data_start.data_arrays[0];
		//prog_u = *parareal_data_start.data_arrays[1];
		//prog_v = *parareal_data_start.data_arrays[2];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if (std::find(this->SL_tsm.begin(), this->SL_tsm.end(), simVars.disc.timestepping_method) != this->SL_tsm.end())
		{
			///PlaneData h_prev = *parareal_data_fine_previous_time_slice.data_arrays[0];
			///PlaneData u_prev = *parareal_data_fine_previous_time_slice.data_arrays[1];
			///PlaneData v_prev = *parareal_data_fine_previous_time_slice.data_arrays[2];

			this->set_previous_solution();
			
			///timeSteppers.master->set_previous_solution(h_prev, u_prev, v_prev);
		}

		while (simVars.timecontrol.current_simulation_time != timeframe_end)
		{
			this->run_timestep(*parareal_data_start, "fine");
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		//// copy to buffers
		//*parareal_data_fine.data_arrays[0] = prog_h_pert;
		//*parareal_data_fine.data_arrays[1] = prog_u;
		//*parareal_data_fine.data_arrays[2] = prog_v;
	};


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_GenericData& get_reference_to_data_timestep_fine()
	{
	};


	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	void run_timestep_coarse()
	{
	};


	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_GenericData& get_reference_to_data_timestep_coarse()
	{
	};

	/**
	 * Compute the error between the fine and coarse timestepping:
	 * Y^E := Y^F - Y^C
	 */
	void compute_difference()
	{
	};

	/**
	 * return the penult solution after the coarse propagation:
	 */
	Parareal_GenericData& get_reference_to_data_timestep_coarse_previous_timestep()
	{
	};

	/**
	 * return the penult solution after the fine propagation:
	 */
	Parareal_GenericData& get_reference_to_data_timestep_fine_previous_timestep()
	{
	};


	/**
	 * Compute the data to be forwarded to the next time step
	 * Y^O := Y^C + Y^E
	 *
	 * Return: If true, the error indicator based on the computed error norm between the
	 * old values and new values
	 */
	double compute_output_data(
			bool i_compute_convergence_test
	)
	{
		double convergence = -1;

		if (!i_compute_convergence_test)
		//if (!i_compute_convergence_test || !output_data_valid)
		{
			*(this->parareal_data_output) = *(this->parareal_data_coarse) - *(this->parareal_data_error);
			//////for (int k = 0; k < 3; k++) {
			//////	*parareal_data_output.data_arrays[k] = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];
			//////	//std::cout << timeframe_end << " " << k << " " << (*parareal_data_output.data_arrays[k] - *parareal_data_fine.data_arrays[k]).reduce_maxAbs() << std::endl;
			//////}

			///////// The following lines are necessary for correctly computing the 1st iteration
			///////// (else, the first time step is not changed from the 0th to the 1st iteration)
			///////// Why?
			///////simVars.timecontrol.current_simulation_time = timeframe_end;
			///////prog_h_pert = *parareal_data_output.data_arrays[0];
			///////prog_u = *parareal_data_output.data_arrays[1];
			///////prog_v = *parareal_data_output.data_arrays[2];

			//output_data_valid = true;
			return convergence;
		}

		Parareal_GenericData* tmp = new t_dataType2<N>;
		*tmp = *(this->parareal_data_coarse) + *(this->parareal_data_error);
		//PlaneData tmp = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

		convergence = (*(this->parareal_data_output)-*tmp).reduce_maxAbs();

		this->parareal_data_output = tmp;
		

		//////////////simVars.timecontrol.current_simulation_time = timeframe_end;
		//////////////prog_h_pert = *parareal_data_output.data_arrays[0];
		//////////////prog_u = *parareal_data_output.data_arrays[1];
		//////////////prog_v = *parareal_data_output.data_arrays[2];

		//////////////if (compute_error_to_analytical_solution)
		//////////////{
		//////////////	if (simVars.misc.compute_errors > 0)
		//////////////	{
		//////////////		compute_errors();
		//////////////		std::cout << "maxabs error compared to analytical solution: " << benchmark.analytical_error_maxabs_h << std::endl;
		//////////////	}
		//////////////}

		//output_data_valid = true;
		return convergence;

	};


	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_GenericData& get_reference_to_output_data()
	{
	};

	void output_data_console(
			int iteration_id,
			int time_slice_id
	)
	{
	};

	void output_data_file(
			int iteration_id,
			int time_slice_id
	)
	{
	};

	void check_for_nan_parareal()
	{
	};

	~Parareal_SimulationInstance_GenericData()
	{
	}
};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_ */
