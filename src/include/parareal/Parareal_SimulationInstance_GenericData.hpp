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

	// Plane
	void setup(SimulationVariables &i_simVars, std::string i_geometry, std::string i_model,
		   PlaneOperators &i_op_plane,
		   t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{
		// call default setup
		this->setup(i_simVars, i_geometry, i_model, i_timeSteppersFine, i_timeSteppersCoarse);
		this->op_plane = i_op_plane;
	}

	// Sphere
	void setup(SimulationVariables &i_simVars, std::string i_geometry, std::string i_model,
		   SphereOperators_SphereData &i_op_sphere, SphereOperators_SphereData &i_op_sphere_nodealiasing,
		   t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{
		// call default setup
		this->setup(i_simVars, i_geometry, i_model, i_timeSteppersFine, i_timeSteppersCoarse);
		this->op_sphere = i_op_sphere;
		this->op_sphere_nodealiasing = i_op_sphere_nodealiasing;
	}


	void setup(SimulationVariables &i_simVars, std::string i_geometry, std::string i_model,
			t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{

		this->timeSteppersFine = i_timeSteppersFine;
		this->timeSteppersCoarse = i_timeSteppersCoarse;

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

	};



public:

	// Functions to exchange data between Parareal_GenericDAta instances and individual data arrays

	void dataArrays_to_GenericData_Scalar(Parareal_GenericData* i_data
						double &u)
	{
		if (this->model == "ode1")
		{
			*(i_data.get_pointer_to_data_Scalar()->simfields[0]) = u;
		}
		else
			SWEETError("Unknown model");
	}

	void GenericData_Scalar_to_DataArrays(Parareal_GenericData* i_data
						double &u)
	{
		if (this->model == "ode1")
		{
			u = *(i_data.get_pointer_to_data_Scalar()->simfields[0]);
		}
		else
			SWEETError("Unknown model");
	}

	void dataArrays_to_GenericData_PlaneData_Spectral(Parareal_GenericData* i_data
						PlaneData &h, PlaneData &u, PlaneData&v)
	{
		if (this->model == "swe")
		{
			*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = h;
			*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = u;
			*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[2]) = v;
		}
		else if (this->model == "burgers")
		{
			*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = u;
			*(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = v;
		}
		else
			SWEETError("Unknown model");
	}

	void GenericData_PlaneData_Spectral_to_dataArrays(Parareal_GenericData* i_pararealData
						PlaneData &h, PlaneData &u, PlaneData&v)
	{
		if (this->model == "swe")
		{
			h = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
			u = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[1]);
			v = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[2]);
		}
		else if (this->model == "burgers")
		{
			u = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
			v = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[1]);
		}
		else
			SWEETError("Unknown model");
	}

	void dataArrays_to_GenericData_SphereData_Spectral(Parareal_GenericData* i_data
						SphereData_Spectral &h, SphereData_Spectral &u, SphereData_Spectral&v)
	{
		if (this->model == "swe")
		{
			*(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[0]) = h;
			*(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[1]) = u;
			*(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[2]) = v;
		}
		else
			SWEETError("Unknown model");
	}

	void GenericData_SphereData_Spectral_to_dataArrays(Parareal_GenericData* i_pararealData
						SphereData_Spectral &h, SphereData_Spectral &u, SphereData_Spectral&v)
	{
		if (this->model == "swe")
		{
			h = *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[0]);
			u = *(i_data.get_pointer_to_data_Data_Spectral()->simfields[1]);
			v = *(i_data.get_pointer_to_data_Data_Spectral()->simfields[2]);
		}
		else
			SWEETError("Unknown model");
	}



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
		if (simVars.parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;
		this->timeframe_start = i_timeframe_start;
		this->timeframe_end = i_timeframe_end;

		// set time to parareal_genericdata instances
		this->parareal_data_start->set_time(i_timeframe_begin);
		this->parareal_data_fine->set_time(i_timeframe_end);
		this->parareal_data_coarse->set_time(i_timeframe_end);
		this->parareal_data_output->set_time(i_timeframe_end);
		this->parareal_data_error->set_time(i_timeframe_end);
		this->parareal_data_coarse_previous_timestep->set_time(i_timeframe_end);
		this->parareal_data_coarse_previous_time_slice->set_time(i_timeframe_end);
		this->parareal_data_fine_previous_timestep->set_time(i_timeframe_end);
		this->parareal_data_fine_previous_time_slice->set_time(i_timeframe_end);
#if SWEET_DEBUG
		this->parareal_data_fine_exact->set_time(i_timeframe_end);
#endif

	};

	/**
	 * Set the initial data at i_timeframe_start
	 */
	void sim_setup_initial_data(
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_setup_initial_data()" << std::endl;

		reset();

		if (this->geometry == "scalar")
		{
			if (this->model == "ode1")
			{
				SWEETError("TODO");
			}
			else
				SWEETError("Unknown model for this geometry");
		}
		else if (this->geometry == "plane")
		{
			PlaneData t0_prog_h_pert(planeDataConfig),
			PlaneData t0_prog_u(planeDataConfig),
			PlaneData t0_prog_v(planeDataConfig),

			if (this->model == "swe")
			{
				SWEPlaneBenchmarksCombined swePlaneBenchmarks;
				swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, simVars, op);
			}
			else if (this->model == "burgers")
			{
				if (simVars.disc.space_grid_use_c_staggering)
				{
					prog_u.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
						}
					);
					prog_v.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
						io_data = 0.0;
						}
					);
				}
				else
				{
					prog_u.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
						}
					);
		
					prog_v.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
						io_data = 0.0;
						}
					);
				}
			}
			else
				SWEETError("Unknown model for this geometry");

			this->dataArrays_to_GenericData_PlaneData_Spectral(this->parareal_data_start, t0_prog_h_pert, t0_prog_u, t0_prog_v)
			this->dataArrays_to_GenericData_PlaneData_Spectral(this->parareal_data_coarse_previous_time_slice, t0_prog_h_pert, t0_prog_u, t0_prog_v)
			this->dataArrays_to_GenericData_PlaneData_Spectral(this->parareal_data_fine_previous_time_slice, t0_prog_h_pert, t0_prog_u, t0_prog_v)

		}
		else if (this->geometry == "sphere")
		{
			PlaneData t0_prog_phi_pert(planeDataConfig),
			PlaneData t0_prog_vrt(planeDataConfig),
			PlaneData t0_prog_div(planeDataConfig),
			if (model == "swe")
			{
				BenchmarksSphereSWE sphereBenchmarks;
				sphereBenchmarks.setup(simVars, op);
				sphereBenchmarks.master->get_initial_state(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
			}
			else
				SWEETError("Unknown model for this geometry");

			this->dataArrays_to_GenericData_SphereData_Spectral(this->parareal_data_start, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div)
			this->dataArrays_to_GenericData_SphereData_Spectral(this->parareal_data_coarse_previous_time_slice, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div)
			this->dataArrays_to_GenericData_SphereData_Spectral(this->parareal_data_fine_previous_time_slice, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div)
		}
		else
			SWEETError("Unknown geometry");
	);

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
	void sim_set_data_fine_previous_time_slice(
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
	void sim_set_data_fine_exact(
			Parareal_GenericData &i_pararealData
	)
	{
		SWEETError("TODO");
	};

	/**
	* Check if solution at time k (end of time slice k-1) is exact (= fine) at iteration k
	*/
	void compare_to_fine_exact(
	)
	{
		SWEETError("TODO");
	};

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
