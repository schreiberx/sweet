/*
 * PararealSimulation.hpp
 *
 *  Created on: 25 Feb 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/BenchmarksSphereSWE.cpp
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_


#include <parareal/Parareal_GenericData.hpp>
#include <parareal/Parareal_GenericData_Scalar.hpp>
#include <parareal/Parareal_GenericData_PlaneData_Spectral.hpp>
#include <parareal/Parareal_GenericData_SphereData_Spectral.hpp>

#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>

#include <sweet/plane/PlaneDataGridMapping.hpp>

#include "../../programs/swe_plane_timeintegrators/SWE_Plane_TimeSteppers.hpp"
#include "../../programs/swe_sphere_timeintegrators/SWE_Sphere_TimeSteppers.hpp"
#include "../../programs/burgers_timeintegrators/Burgers_Plane_TimeSteppers.hpp"

#include "../../programs/swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"
#include "../../programs/swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
//#include "../../programs/swe_sphere_benchmarks/BenchmarksSphereSWE.cpp"

/**
 * Interface descriptions which are required
 * to run Parareal Simulations in SWEET.
 *
 * Note that these interface descriptions have
 * to be implemented by the simulation!
 *
 * These interfaces were ported from the Python implementation.
 */
template <class t_tsmType, int N>
class Parareal_SimulationInstance_GenericData
{
public:

	// Simulation variables
	SimulationVariables* simVars;

	// Grid Mapping (staggered grid)
	PlaneDataGridMapping gridMapping;

	// Operators
	PlaneOperators* op_plane;
	SphereOperators_SphereData* op_sphere;
	SphereOperators_SphereData* op_sphere_nodealiasing;

	// Data config
	PlaneDataConfig* planeDataConfig;
	SphereData_Config* sphereDataConfig;

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
	Parareal_GenericData* parareal_data_ref_exact = nullptr;
//#if SWEET_DEBUG
	Parareal_GenericData* parareal_data_fine_exact = nullptr;
//#endif


	// Fine and coarse timesteppers
	t_tsmType* timeSteppersFine = nullptr;
	t_tsmType* timeSteppersCoarse = nullptr;

	// list of SL schemes
	std::vector<std::string> SL_tsm = {};


	bool compute_normal_modes = false;

public:

	Parareal_SimulationInstance_GenericData()
	{
	};

	// Plane
	void setup(SimulationVariables* i_simVars, PlaneDataConfig* i_planeDataConfig,
		   std::string i_geometry, std::string i_model,
		   PlaneOperators* i_op_plane,
		   t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{
		this->planeDataConfig = i_planeDataConfig;
		this->op_plane = i_op_plane;
		this->setup(i_simVars, i_geometry, i_model, i_timeSteppersFine, i_timeSteppersCoarse);

		this->SL_tsm = { "l_cn_na_sl_nd_settls",
				 "l_rexi_na_sl_nd_etdrk",
				 "l_rexi_na_sl_nd_settls"
				};

		if (simVars->benchmark.benchmark_name == "normalmodes" )
			this->compute_normal_modes = true;

	}

	// Sphere
	void setup(SimulationVariables* i_simVars, SphereData_Config* i_sphereDataConfig,
		   std::string i_geometry, std::string i_model,
		   SphereOperators_SphereData* i_op_sphere, SphereOperators_SphereData* i_op_sphere_nodealiasing,
		   t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{
		this->sphereDataConfig = i_sphereDataConfig;
		this->op_sphere = i_op_sphere;
		this->op_sphere_nodealiasing = i_op_sphere_nodealiasing;
		this->setup(i_simVars, i_geometry, i_model, i_timeSteppersFine, i_timeSteppersCoarse);

		this->SL_tsm = { "lg_exp_na_sl_lc_nr_etd_uv",
				 "l_irk_na_sl_nr_settls_uv_only",
				 "l_irk_na_sl_nr_settls_vd_only",
				 "l_irk_na_sl_settls_uv_only",
				 "l_irk_na_sl_settls_vd_only",
				 "ln_sl_exp_settls_uv",
				 "ln_sl_exp_settls_vd"
				};
	}


	void setup(SimulationVariables* i_simVars, std::string i_geometry, std::string i_model,
			t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{

		this->simVars = i_simVars;

		this->geometry = i_geometry;
		this->model = i_model;

		this->timeSteppersFine = i_timeSteppersFine;
		this->timeSteppersCoarse = i_timeSteppersCoarse;

		this->parareal_data_start                      = this->create_new_data_container();
		this->parareal_data_fine                       = this->create_new_data_container();
		this->parareal_data_coarse                     = this->create_new_data_container();
		this->parareal_data_output                     = this->create_new_data_container();
		this->parareal_data_error                      = this->create_new_data_container();
		this->parareal_data_coarse_previous_timestep   = this->create_new_data_container();
		this->parareal_data_coarse_previous_time_slice = this->create_new_data_container();
		this->parareal_data_fine_previous_timestep     = this->create_new_data_container();
		this->parareal_data_fine_previous_time_slice   = this->create_new_data_container();
		this->parareal_data_ref_exact   = this->create_new_data_container();
//#if SWEET_DEBUG
		this->parareal_data_fine_exact                 = this->create_new_data_container();
//#endif
	};


public:

	Parareal_GenericData* create_new_data_container()
	{
		if (this->geometry == "scalar") 
		{
			Parareal_GenericData_Scalar<N>* out = new Parareal_GenericData_Scalar<N>;
			out->allocate_data();
			return out;
		}
		else if (this->geometry == "plane")
		{
			Parareal_GenericData_PlaneData_Spectral<N>* out = new Parareal_GenericData_PlaneData_Spectral<N>;
			out->setup_data_config(this->planeDataConfig);
			out->allocate_data();
			return out;
		}
		else if (this->geometry == "sphere")
		{
			Parareal_GenericData_SphereData_Spectral<N>* out = new Parareal_GenericData_SphereData_Spectral<N>;
			out->setup_data_config(this->sphereDataConfig);
			out->allocate_data();
			return out;
		}
		else
			SWEETError("Unknown geometry");
	}



public:

	// Functions to exchange data between Parareal_GenericData instances and individual data arrays

	void dataArrays_to_GenericData_Scalar(Parareal_GenericData* i_data,
						double &u)
	{
		if (this->model == "ode1")
		{
			i_data->get_pointer_to_data_Scalar()->simfields[0] = u;
		}
		else
			SWEETError("Unknown model");
	}

	void GenericData_Scalar_to_DataArrays(Parareal_GenericData* i_data,
						double &u)
	{
		if (this->model == "ode1")
		{
			u = i_data->get_pointer_to_data_Scalar()->simfields[0];
		}
		else
			SWEETError("Unknown model");
	}

	void dataArrays_to_GenericData_PlaneData_Spectral(Parareal_GenericData* i_data,
						PlaneData_Spectral &h, PlaneData_Spectral &u, PlaneData_Spectral &v)
	{
		if (this->model == "swe")
		{
			*(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = h;
			*(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = u;
			*(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2]) = v;
		}
		else if (this->model == "burgers")
		{
			*(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = u;
			*(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = v;
		}
		else
			SWEETError("Unknown model");
	}

	void GenericData_PlaneData_Spectral_to_dataArrays(Parareal_GenericData* i_data,
						PlaneData_Spectral &h, PlaneData_Spectral &u, PlaneData_Spectral &v)
	{
		if (this->model == "swe")
		{
			h = *(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
			u = *(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]);
			v = *(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2]);
		}
		else if (this->model == "burgers")
		{
			u = *(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
			v = *(i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]);
		}
		else
			SWEETError("Unknown model");
	}

	void dataArrays_to_GenericData_SphereData_Spectral(Parareal_GenericData* i_data,
						SphereData_Spectral &phi, SphereData_Spectral &vrt, SphereData_Spectral& div)
	{
		if (this->model == "swe")
		{
			*(i_data->get_pointer_to_data_SphereData_Spectral()->simfields[0]) = phi;
			*(i_data->get_pointer_to_data_SphereData_Spectral()->simfields[1]) = vrt;
			*(i_data->get_pointer_to_data_SphereData_Spectral()->simfields[2]) = div;
		}
		else
			SWEETError("Unknown model");
	}

	void GenericData_SphereData_Spectral_to_dataArrays(Parareal_GenericData* i_data,
						SphereData_Spectral &phi, SphereData_Spectral &vrt, SphereData_Spectral& div)
	{
		if (this->model == "swe")
		{
			phi = *(i_data->get_pointer_to_data_SphereData_Spectral()->simfields[0]);
			vrt = *(i_data->get_pointer_to_data_SphereData_Spectral()->simfields[1]);
			div = *(i_data->get_pointer_to_data_SphereData_Spectral()->simfields[2]);
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
		double mod_coarse = fmod(time_slice_size, simVars->parareal.coarse_timestep_size);
		double mod_fine = fmod(time_slice_size, simVars->timecontrol.current_timestep_size);
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
		if (simVars->parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;
		this->timeframe_start = i_timeframe_start;
		this->timeframe_end = i_timeframe_end;

		// set time to parareal_genericdata instances
		this->parareal_data_start->set_time(i_timeframe_end);
		this->parareal_data_fine->set_time(i_timeframe_end);
		this->parareal_data_coarse->set_time(i_timeframe_end);
		this->parareal_data_output->set_time(i_timeframe_end);
		this->parareal_data_error->set_time(i_timeframe_end);
		this->parareal_data_coarse_previous_timestep->set_time(i_timeframe_end);
		this->parareal_data_coarse_previous_time_slice->set_time(i_timeframe_end);
		this->parareal_data_fine_previous_timestep->set_time(i_timeframe_end);
		this->parareal_data_fine_previous_time_slice->set_time(i_timeframe_end);
		this->parareal_data_ref_exact->set_time(i_timeframe_end);
		this->parareal_data_fine_exact->set_time(i_timeframe_end);

	};

	/**
	 * Set the initial data at i_timeframe_start
	 */
	void sim_setup_initial_data()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_setup_initial_data()" << std::endl;

		///reset();

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
			PlaneData_Spectral t0_prog_h_pert(planeDataConfig);
			PlaneData_Spectral t0_prog_u(planeDataConfig);
			PlaneData_Spectral t0_prog_v(planeDataConfig);

			if (this->model == "swe")
			{
				SWEPlaneBenchmarksCombined swePlaneBenchmarks;
				swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, *simVars, *op_plane);
			}
			else if (this->model == "burgers")
			{
				PlaneData_Physical t0_prog_u_phys(t0_prog_u.planeDataConfig);
				PlaneData_Physical t0_prog_v_phys(t0_prog_v.planeDataConfig);
				if (simVars->disc.space_grid_use_c_staggering)
				{
					t0_prog_u_phys.physical_update_lambda_array_indices(
								[&](int i, int j, double &io_data)
						{
							double x = (((double)i)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);
					t0_prog_v_phys.physical_update_lambda_array_indices(
								[&](int i, int j, double &io_data)
						{
						io_data = 0.0;
						}
					);
				}
				else
				{
					t0_prog_u_phys.physical_update_lambda_array_indices(
								[&](int i, int j, double &io_data)
						{
							double x = (((double)i+0.5)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);
		
					t0_prog_v_phys.physical_update_lambda_array_indices(
								[&](int i, int j, double &io_data)
						{
						io_data = 0.0;
						}
					);
				}
				t0_prog_u.loadPlaneDataPhysical(t0_prog_u_phys);
				t0_prog_v.loadPlaneDataPhysical(t0_prog_v_phys);
			}
			else
				SWEETError("Unknown model for this geometry");

			this->dataArrays_to_GenericData_PlaneData_Spectral(this->parareal_data_start, t0_prog_h_pert, t0_prog_u, t0_prog_v);
			this->dataArrays_to_GenericData_PlaneData_Spectral(this->parareal_data_coarse_previous_time_slice, t0_prog_h_pert, t0_prog_u, t0_prog_v);
			this->dataArrays_to_GenericData_PlaneData_Spectral(this->parareal_data_fine_previous_time_slice, t0_prog_h_pert, t0_prog_u, t0_prog_v);
		}
		else if (this->geometry == "sphere")
		{
			SphereData_Spectral t0_prog_phi_pert(sphereDataConfig);
			SphereData_Spectral t0_prog_vrt(sphereDataConfig);
			SphereData_Spectral t0_prog_div(sphereDataConfig);
			if (model == "swe")
			{
				BenchmarksSphereSWE sphereBenchmarks;
				sphereBenchmarks.setup(*simVars, *op_sphere);
				sphereBenchmarks.master->get_initial_state(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
			}
			else
				SWEETError("Unknown model for this geometry");

			this->dataArrays_to_GenericData_SphereData_Spectral(this->parareal_data_start, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
			this->dataArrays_to_GenericData_SphereData_Spectral(this->parareal_data_coarse_previous_time_slice, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
			this->dataArrays_to_GenericData_SphereData_Spectral(this->parareal_data_fine_previous_time_slice, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
		}
		else
			SWEETError("Unknown geometry");
	};


	/**
	 * Set simulation data to data given in i_sim_data.
	 * This can be data which is computed by another simulation.
	 * Y^S := i_sim_data
	 */
	void sim_set_data(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
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
		if (simVars->parareal.verbosity > 2)
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
		if (simVars->parareal.verbosity > 2)
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
		//SWEETError("TODO");
	};

	/**
	* Check if solution at time k (end of time slice k-1) is exact (= fine) at iteration k
	*/
	void compare_to_fine_exact(
	)
	{
		//SWEETError("TODO");
	};

#endif

	void set_previous_solution(
				std::string tsm_level
				)
	{

		if (tsm_level == "fine")
			timeSteppersFine->master->set_previous_solution(this->parareal_data_fine_previous_time_slice);
		else if (tsm_level == "coarse")
			timeSteppersCoarse->master->set_previous_solution(this->parareal_data_fine_previous_time_slice);
		else
			SWEETError("Wrong tsm_level (should be 'fine' or 'coarse')");

	};

	/**
	 * Set the MPI communicator to use for simulation purpose
	 * (TODO: not yet implemented since our parallelization-in-space
	 * is done only via OpenMP)
	 */
	void sim_set_mpi_comm(
			int i_mpi_comm
	)
	{
	};


	void run_timestep(
			Parareal_GenericData* i_data,
			std::string tsm_level
	)
	{

		if (tsm_level == "fine")
			timeSteppersFine->master->run_timestep(
						i_data,
						simVars->timecontrol.current_timestep_size,
						simVars->timecontrol.current_simulation_time
					);
		else if (tsm_level == "coarse")
			timeSteppersCoarse->master->run_timestep(
						i_data,
						simVars->timecontrol.current_timestep_size,
						simVars->timecontrol.current_simulation_time
					);
		else
			SWEETError("Wrong tsm_level (should be 'fine' or 'coarse')");

	}



	/**
	 * compute solution on time slice with fine timestep:
	 * Y^F := F(Y^S)
	 */
	void run_timestep_fine()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "run_timestep_fine()" << std::endl;

		// reset simulation time
		simVars->timecontrol.current_simulation_time = timeframe_start;
		simVars->timecontrol.max_simulation_time = timeframe_end;
		simVars->timecontrol.current_timestep_nr = 0;

		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if (std::find(this->SL_tsm.begin(), this->SL_tsm.end(), simVars->disc.timestepping_method) != this->SL_tsm.end())
			this->set_previous_solution("fine");

		*(this->parareal_data_fine) = *(this->parareal_data_start);

		while (simVars->timecontrol.current_simulation_time != timeframe_end)
		{
			// store previous time step
			// to be used as n-1 in SL in the next time slice
			*(this->parareal_data_fine_previous_timestep) = *(this->parareal_data_fine);

			this->run_timestep(this->parareal_data_fine, "fine");

			simVars->timecontrol.current_simulation_time += simVars->timecontrol.current_timestep_size;
			assert(simVars->timecontrol.current_simulation_time <= timeframe_end);
		}
	};


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_GenericData& get_reference_to_data_timestep_fine()
	{
		return *(this->parareal_data_fine);
	};


	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	void run_timestep_coarse()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "run_timestep_coarse()" << std::endl;

		// reset simulation time
		simVars->timecontrol.current_simulation_time = timeframe_start;
		simVars->timecontrol.max_simulation_time = timeframe_end;
		simVars->timecontrol.current_timestep_nr = 0;

		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if (std::find(this->SL_tsm.begin(), this->SL_tsm.end(), simVars->parareal.coarse_timestepping_method) != this->SL_tsm.end())
			this->set_previous_solution("coarse");

		*(this->parareal_data_coarse) = *(this->parareal_data_start);

		while (simVars->timecontrol.current_simulation_time != timeframe_end)
		{
			// store previous time step
			// to be used as n-1 in SL in the next time slice
			*(this->parareal_data_coarse_previous_timestep) = *(this->parareal_data_coarse);

			this->run_timestep(this->parareal_data_coarse, "coarse");
			simVars->timecontrol.current_simulation_time += simVars->parareal.coarse_timestep_size;
			assert(simVars->timecontrol.current_simulation_time <= timeframe_end);
		}
	};


	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_GenericData& get_reference_to_data_timestep_coarse()
	{
		return *(this->parareal_data_coarse);
	};

	/**
	 * Compute the error between the fine and coarse timestepping:
	 * Y^E := Y^F - Y^C
	 */
	void compute_difference()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "compute_difference()" << std::endl;

		*(this->parareal_data_error) = *(this->parareal_data_fine);
		*(this->parareal_data_error) -= *(this->parareal_data_coarse);
	};

	/**
	 * return the penult solution after the coarse propagation:
	 */
	Parareal_GenericData& get_reference_to_data_timestep_coarse_previous_timestep()
	{
		return *(this->parareal_data_coarse_previous_timestep);
	};

	/**
	 * return the penult solution after the fine propagation:
	 */
	Parareal_GenericData& get_reference_to_data_timestep_fine_previous_timestep()
	{
		return *(this->parareal_data_fine_previous_timestep);
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
			*(this->parareal_data_output) = *(this->parareal_data_coarse);
                        *(this->parareal_data_output) -= *(this->parareal_data_error);

			//output_data_valid = true;
			return convergence;
		}

		// compute output data
		Parareal_GenericData* tmp = this->create_new_data_container();
		*tmp = *(this->parareal_data_coarse);
                *tmp += *(this->parareal_data_error);

		// compute difference w.r.t. previous output data
		Parareal_GenericData* tmp2 = this->create_new_data_container();
		*tmp2 = *(this->parareal_data_output);
		*tmp2 -= *tmp;
		convergence = tmp2->reduce_maxAbs();

		// store output data
		*(this->parareal_data_output) = *tmp;

		delete tmp;
		delete tmp2;

		//output_data_valid = true;
		return convergence;

	};


	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_GenericData& get_reference_to_output_data()
	{
		return *(this->parareal_data_output);
	};

	void output_data_console(
			int iteration_id,
			int time_slice_id
	)
	{
	};

	void output_data_file(
			int iteration_id,
			int time_slice_id,
			bool output_initial_data = false
	)
	{
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

				PlaneData_Spectral h_out(this->planeDataConfig);
				PlaneData_Spectral u_out(this->planeDataConfig);
				PlaneData_Spectral v_out(this->planeDataConfig);
				if (output_initial_data)
					this->GenericData_PlaneData_Spectral_to_dataArrays(this->parareal_data_start, h_out, u_out, v_out);
				else
					this->GenericData_PlaneData_Spectral_to_dataArrays(this->parareal_data_output, h_out, u_out, v_out);

				PlaneData_Physical h_out_phys = h_out.toPhys();
				PlaneData_Physical u_out_phys = u_out.toPhys();
				PlaneData_Physical v_out_phys = v_out.toPhys();

				// Save .vtk files for visualizing in paraview
				std::ostringstream ss2;
				if (output_initial_data)
					ss2 << "output_slice" << time_slice_id - 1 << "_iter" << iteration_id << ".vtk";
				else
					ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
				std::string filename2 = ss2.str();
				h_out_phys.file_physical_saveData_vtk(filename2.c_str(), filename2.c_str());

				/*
				 * File output
				 *
				 * We write everything in non-staggered output
				 */
				// For output, variables need to be on unstaggered A-grid
				PlaneData_Physical t_h(planeDataConfig);
				PlaneData_Physical t_u(planeDataConfig);
				PlaneData_Physical t_v(planeDataConfig);

				if (simVars->disc.space_grid_use_c_staggering) // Remap in case of C-grid
				{
					t_h = h_out_phys;
					gridMapping.mapCtoA_u(u_out_phys, t_u);
					gridMapping.mapCtoA_v(v_out_phys, t_v);
				}
				else
				{
					t_h = h_out_phys;
					t_u = u_out_phys;
					t_v = v_out_phys;
				}

				// Dump  data in csv, if output filename is not empty
				if (simVars->iodata.output_file_name.size() > 0)
				{
					std::string output_filenames = "";

					output_filenames = write_file_parareal_plane(t_h, "prog_h_pert", iteration_id, output_initial_data);
					output_filenames += ";" + write_file_parareal_plane(t_u, "prog_u", iteration_id, output_initial_data);
					output_filenames += ";" + write_file_parareal_plane(t_v, "prog_v", iteration_id, output_initial_data);

					output_filenames += ";" + write_file_parareal_plane(op_plane->ke(t_u,t_v).toPhys(),"diag_ke", iteration_id, output_initial_data);

					output_filenames += ";" + write_file_spec_parareal_plane(op_plane->ke(t_u,t_v).toPhys(),"diag_ke_spec", iteration_id, output_initial_data);

					output_filenames += ";" + write_file_parareal_plane(op_plane->vort(t_u, t_v).toPhys(), "diag_vort", iteration_id, output_initial_data);
					output_filenames += ";" + write_file_parareal_plane(op_plane->div(t_u, t_v).toPhys(), "diag_div", iteration_id, output_initial_data);

					if(this->compute_normal_modes){
						SWEETError("TODO");
						///output_filenames += ";" + write_file_spec_parareal(normalmodes.geo, "nm_geo", iteration_id, output_initial_data);
						///output_filenames += ";" + write_file_spec_parareal(normalmodes.igwest, "nm_igwest", iteration_id, output_initial_data);
						///output_filenames += ";" + write_file_spec_parareal(normalmodes.igeast, "nm_igeast", iteration_id, output_initial_data);
					}
					
				}
			}
		}
		else if (this->geometry == "sphere")
		{

			if (this->model == "swe")
			{
				SphereData_Spectral phi_out(this->sphereDataConfig);
				SphereData_Spectral vrt_out(this->sphereDataConfig);
				SphereData_Spectral div_out(this->sphereDataConfig);
				if (output_initial_data)
					this->GenericData_SphereData_Spectral_to_dataArrays(this->parareal_data_start, phi_out, vrt_out, div_out);
				else
					this->GenericData_SphereData_Spectral_to_dataArrays(this->parareal_data_output, phi_out, vrt_out, div_out);
	
				SphereData_Physical phi_out_phys = phi_out.toPhys();
				SphereData_Physical vrt_out_phys = vrt_out.toPhys();
				SphereData_Physical div_out_phys = div_out.toPhys();
	
				///////////////// Save .vtk files for visualizing in paraview
				///////////////std::ostringstream ss2;
				///////////////if (output_initial_data)
				///////////////	ss2 << "output_slice" << time_slice_id - 1 << "_iter" << iteration_id << ".vtk";
				///////////////else
				///////////////	ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
				///////////////std::string filename2 = ss2.str();
				///////////////phi_out_phys.file_physical_saveData_vtk(filename2.c_str(), filename2.c_str());


				/*
				 * File output
				 *
				 * We write everything in non-staggered output
				 */
				// Dump  data in csv, if output filename is not empty
				if (simVars->iodata.output_file_name.size() > 0)
				{
					if (simVars->iodata.output_file_mode == "csv")
					{
						std::string output_filename;
		
						SphereData_Spectral h = phi_out_phys*(1.0/simVars->sim.gravitation);
						h += simVars->sim.h0;
		
						output_filename = write_file_csv_parareal_sphere(h, "prog_h", iteration_id, output_initial_data);
						std::cout << " + " << output_filename << " (min: " << h.toPhys().physical_reduce_min() << ", max: " << h.toPhys().physical_reduce_max() << ")" << std::endl;
		
						output_filename = write_file_csv_parareal_sphere(phi_out, "prog_phi_pert", iteration_id, output_initial_data);
						std::cout << " + " << output_filename << " (min: " << phi_out_phys.physical_reduce_min() << ", max: " << phi_out_phys.physical_reduce_max() << ")" << std::endl;
		
						SphereData_Physical u(sphereDataConfig);
						SphereData_Physical v(sphereDataConfig);
		
						op_sphere->vrtdiv_to_uv(vrt_out_phys, div_out_phys, u, v);
		
						output_filename = write_file_csv_parareal_sphere(u, "prog_u", iteration_id, output_initial_data);
						std::cout << " + " << output_filename << std::endl;
		
						output_filename = write_file_csv_parareal_sphere(v, "prog_v", iteration_id, output_initial_data);
						std::cout << " + " << output_filename << std::endl;
		
						output_filename = write_file_csv_parareal_sphere(vrt_out, "prog_vrt", iteration_id, output_initial_data);
						std::cout << " + " << output_filename << std::endl;
		
						output_filename = write_file_csv_parareal_sphere(div_out, "prog_div", iteration_id, output_initial_data);
						std::cout << " + " << output_filename << std::endl;
		
						SphereData_Spectral potvrt = (phi_out/simVars->sim.gravitation)*vrt_out;
		
						output_filename = write_file_csv_parareal_sphere(potvrt, "prog_potvrt", iteration_id, output_initial_data);
						std::cout << " + " << output_filename << std::endl;
					}
					else if (simVars->iodata.output_file_mode == "bin")
					{
						std::string output_filename;
		
						{
							output_filename = write_file_bin_parareal_sphere(phi_out, "prog_phi_pert", iteration_id);
							SphereData_Physical prog_phys = phi_out.toPhys();
		
							std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
						}
		
						{
							output_filename = write_file_bin_parareal_sphere(vrt_out, "prog_vrt", iteration_id);
							SphereData_Physical prog_phys = vrt_out.toPhys();
		
							std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
						}
		
						{
							output_filename = write_file_bin_parareal_sphere(div_out, "prog_div", iteration_id);
							SphereData_Physical prog_phys = div_out.toPhys();
		
							std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
						}
					}
					
				}

			}

		}
		else
			SWEETError("Unknown geometry");
	};


	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_csv_parareal_sphere(
			const SphereData_Spectral &i_sphereData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			bool output_initial_data = false,
			bool i_phi_shifted = false
		)
	{
		char buffer[1024];

		// create copy
		SphereData_Physical sphereData = i_sphereData.toPhys();

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, timeframe_end * simVars->iodata.output_time_scale, iteration_id);

		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;

	}

	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_bin_parareal_sphere(
			const SphereData_Spectral &i_sphereData,
			const char* i_name,
			int iteration_id
	)
	{
		char buffer[1024];

		SphereData_Spectral sphereData(i_sphereData);
		//const char* filename_template = simVars.iodata.output_file_name.c_str();
		const char* filename_template = "output_%s_t%020.8f_iter%03d.sweet";
		sprintf(buffer, filename_template, i_name, timeframe_end * simVars->iodata.output_time_scale, iteration_id);
		//sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);
		sphereData.file_write_binary_spectral(buffer);

		return buffer;
	}



	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_parareal_plane(
			const PlaneData_Physical &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			bool output_initial_data = false
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, timeframe_start, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, timeframe_end, iteration_id);
		i_planeData.file_physical_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */

	std::string write_file_spec_parareal_plane(
			const PlaneData_Spectral &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			bool output_initial_data = false
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, timeframe_start, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, timeframe_end, iteration_id);
		i_planeData.file_spectral_abs_saveData_ascii(buffer);
		return buffer;
	}


	/**
	 * Compute and store parareal errors during simulation
	 */
	void store_parareal_error(
			int iteration_id,
			int time_slice_id,
			std::string path_ref,
			std::string base_solution,	// "ref" or "fine"
			int i_precision = 16
	)
	{

		Parareal_GenericData* parareal_data_ref;
		if (base_solution == "ref")
			parareal_data_ref = this->parareal_data_ref_exact;
		else if (base_solution == "fine")
			parareal_data_ref = this->parareal_data_fine_exact;
		else
			SWEETError("Wrong base solution for computing parareal errors.");

		// READ REF DATA
		//Parareal_GenericData& data;
		//if (output_initial_data)
		//	data = (Parareal_GenericData&)this->parareal_data_start;
		//else
		//	data = (Parareal_GenericData&)this->parareal_data_output;
		//Parareal_GenericData& data = (Parareal_GenericData&)this->parareal_data_output;

		if (this->geometry == "scalar")
		{
			SWEETError("TODO");
		}
		else if (this->geometry == "plane")
		{
			PlaneData_Spectral ref_data[] = { PlaneData_Spectral(this->planeDataConfig),
					                  PlaneData_Spectral(this->planeDataConfig),
					                  PlaneData_Spectral(this->planeDataConfig)};

			for (int ivar = 0; ivar < 3; ivar++)
			{
				std::string i_name;
				if (ivar == 0)
					i_name = "prog_h_pert";
				else if (ivar == 1)
					i_name = "prog_u";
				else if (ivar == 2)
					i_name = "prog_v";

				if (iteration_id == 0)
				{
					// load ref file
					char buffer[1024];
					const char* filename_template = simVars->iodata.output_file_name.c_str();
					sprintf(buffer, filename_template, i_name.c_str(), timeframe_end);
					std::string buffer2 = path_ref + "/" + std::string(buffer);
                                        PlaneData_Physical tmp(this->planeDataConfig);
					tmp.file_physical_loadRefData_Parareal(buffer2.c_str());
					ref_data[ivar].loadPlaneDataPhysical(tmp);

					// If necessary, interpolate to coarsest spatial grid
					if (	this->planeDataConfig->physical_res[0] != ref_data[ivar].planeDataConfig->physical_res[0] ||
						this->planeDataConfig->physical_res[1] != ref_data[ivar].planeDataConfig->physical_res[1]
					)
					{
						//TODO
	///					/*
	///					 * setup sampler
	///					 */
	///					PlaneDataSampler sampler2D;
	///					sampler2D.setup(simVars.sim.plane_domain_size, planeDataConfig);
	///		
	///		
	///						/*
	///						 * sample with BiLinear interpolation
	///						 */
	///						PlaneData prog_h3_bilinear(planeDataConfig3);
	///		
	///						sampler2D.bilinear_scalar(
	///								prog_h_pert,	///< input scalar field
	///								Convert_PlaneData_To_ScalarDataArray::physical_convert(px),
	///								Convert_PlaneData_To_ScalarDataArray::physical_convert(py),
	///								prog_h3_bilinear
	///						);
					}
					this->dataArrays_to_GenericData_PlaneData_Spectral(parareal_data_ref, ref_data[0], ref_data[1], ref_data[2]);
				}
			}
			
		}
		else if (this->geometry == "sphere")
		{
		}

		// COMPUTE AND STORE ERRORS
		for (int ivar = 0; ivar < 3; ivar++)
		{

			int resx_data;
			int resy_data;
			double err_L1;
			double err_L2;
			double err_Linf;
			std::string i_name;

			if (this->geometry == "scalar")
			{
				SWEETError("TODO");
			}
			else if (this->geometry == "plane")
			{

				if (ivar == 0)
					i_name = "prog_h_pert";
				else if (ivar == 1)
					i_name = "prog_u";
				else if (ivar == 2)
					i_name = "prog_v";

				resx_data = this->planeDataConfig->physical_res[0];
				resy_data = this->planeDataConfig->physical_res[1];

				PlaneData_Physical diff = this->parareal_data_output->get_pointer_to_data_PlaneData_Spectral()->simfields[ivar]->toPhys() -
                                                          parareal_data_ref->get_pointer_to_data_PlaneData_Spectral()->simfields[ivar]->toPhys();
				err_L1 = diff.physical_reduce_norm1() / (resx_data * resy_data);
				err_L2 = diff.physical_reduce_norm2() / std::sqrt(resx_data * resy_data);
				err_Linf = diff.physical_reduce_max_abs();

			}
			else if (this->geometry == "sphere")
			{
			}

			// save errors in file
			char buffer_out[1024];

			const char* filename_template_out = "parareal_error_%s_%s_t%020.8f_iter%03d.csv";
			sprintf(buffer_out, filename_template_out, base_solution.c_str(), i_name.c_str(), timeframe_end, iteration_id);

			std::ofstream file(buffer_out, std::ios_base::trunc);
			file << std::setprecision(i_precision);

			file << "#BASESOLUTION " << base_solution << " " << path_ref << std::endl;
			file << "#VAR " << i_name << std::endl;
			file << "#ITERATION " << iteration_id << std::endl;
			file << "#TIMESLICE " << time_slice_id << std::endl;
			file << "#TIMEFRAMEEND " << timeframe_end << std::endl;
			file << "errL1 " << err_L1 << std::endl;
			file << "errL2 " << err_L2 << std::endl;
			file << "errLinf " << err_Linf << std::endl;

			file.close();
		}
		parareal_data_ref = nullptr;
	}

	void check_for_nan_parareal()
	{
		if (this->parareal_data_output->check_for_nan())
			SWEETError("Instability detected in parareal!");
	};

	~Parareal_SimulationInstance_GenericData()
	{
	}
};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_GENERICDATA_HPP_ */
