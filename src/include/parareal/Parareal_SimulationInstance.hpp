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

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_HPP_



#include <common_pint/PInT_Common.hpp>

#include <sweet/SimulationVariables.hpp>

#include <parareal/Parareal_GenericData.hpp>

#if SWEET_PARAREAL_SCALAR
#include <parareal/Parareal_GenericData_Scalar.hpp>

#elif SWEET_PARAREAL_PLANE
#include <parareal/Parareal_GenericData_PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataGridMapping.hpp>
#include "../../programs/swe_plane_timeintegrators/SWE_Plane_TimeSteppers.hpp"
#include "../../programs/burgers_timeintegrators/Burgers_Plane_TimeSteppers.hpp"
#include "../../programs/swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"

#elif SWEET_PARAREAL_SPHERE
#include <parareal/Parareal_GenericData_SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include "../../programs/swe_sphere_timeintegrators/SWE_Sphere_TimeSteppers.hpp"
#include "../../programs/swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
#endif
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
class Parareal_SimulationInstance
			: public PInT_Common
{
public:

	// Simulation variables coarse level
	SimulationVariables* simVars_coarse = nullptr; // coarse level (dt, tsm, tso provided by specific parareal parameters)

	// Time slice
	double timeframe_start;
	double timeframe_end;
	int nb_timesteps_fine;
	int nb_timesteps_coarse;
	double dt_fine;
	double dt_coarse;

	// Data containers
	Parareal_GenericData* parareal_data_start = nullptr;
	Parareal_GenericData* parareal_data_fine = nullptr;
	Parareal_GenericData* parareal_data_coarse = nullptr; // defined on fine spatial mesh; used almost everywhere
	Parareal_GenericData* parareal_data_coarse_coarse_mesh = nullptr; // defined on coarse spatial mesh; used only on run_timestep_coarse
	Parareal_GenericData* parareal_data_output = nullptr;
	Parareal_GenericData* parareal_data_error = nullptr;
	Parareal_GenericData* parareal_data_coarse_previous_timestep = nullptr;
	Parareal_GenericData* parareal_data_coarse_previous_time_slice = nullptr;
	Parareal_GenericData* parareal_data_fine_previous_timestep = nullptr;
	Parareal_GenericData* parareal_data_fine_previous_time_slice = nullptr;
	Parareal_GenericData* parareal_data_ref_exact = nullptr;
	Parareal_GenericData* parareal_data_fine_exact = nullptr;
#if SWEET_DEBUG
	Parareal_GenericData* parareal_data_fine_exact_debug = nullptr;
#endif

	Parareal_GenericData* parareal_data_debug = nullptr;
	bool debug_contains_data = false;

	// Fine and coarse timesteppers
	t_tsmType* timeSteppersFine = nullptr;
	t_tsmType* timeSteppersCoarse = nullptr;

	///// list of SL schemes
	///std::vector<std::string> SL_tsm = {};


	bool compute_normal_modes = false;


	// For Burgers
	class BenchmarkErrors
	{
	public:
		// Max difference to initial conditions
		double benchmark_diff_u;
		double benchmark_diff_v;

		// Error measures L2 norm
		double benchmark_analytical_error_rms_u;
		double benchmark_analytical_error_rms_v;

		// Error measures max norm
		double benchmark_analytical_error_maxabs_u;
		double benchmark_analytical_error_maxabs_v;
	};

	BenchmarkErrors benchmark;



public:

	Parareal_SimulationInstance()
	{
	};

#if SWEET_PARAREAL_PLANE
	// Plane
	void setup(	SimulationVariables* i_simVars,
			std::vector<PlaneDataConfig*> i_planeDataConfig,
			std::vector<PlaneOperators*> i_op_plane,
			t_tsmType* i_timeSteppersFine,
			t_tsmType* i_timeSteppersCoarse)
	{
		this->planeDataConfig = i_planeDataConfig;
		this->op_plane = i_op_plane;
		this->setup(i_simVars, 
				i_timeSteppersFine, i_timeSteppersCoarse);

		// IMPORTANT: setup initial conditions (inside this->setup) before setting up timesteppers
		// because simulation parameters may change
		this->timeSteppersFine->setup(
				this->simVars->disc.timestepping_method,
				this->simVars->disc.timestepping_order,
				this->simVars->disc.timestepping_order2,
				*this->op_plane[0],
				*this->simVars
			);

		this->timeSteppersCoarse->setup(
				this->simVars_coarse->disc.timestepping_method,
				this->simVars_coarse->disc.timestepping_order,
				this->simVars_coarse->disc.timestepping_order2,
				*this->op_plane[1],
				*this->simVars_coarse
			);

		PInT_Common::setup();

		if (simVars->benchmark.benchmark_name == "normalmodes" )
			this->compute_normal_modes = true;

	}

#elif SWEET_PARAREAL_SPHERE
	// Sphere
	void setup(	SimulationVariables* i_simVars,
			std::vector<SphereData_Config*> i_sphereDataConfig,
			std::vector<SphereOperators_SphereData*> i_op_sphere,
			std::vector<SphereOperators_SphereData*> i_op_sphere_nodealiasing,
			t_tsmType* i_timeSteppersFine,
			t_tsmType* i_timeSteppersCoarse)
	{
		this->sphereDataConfig = i_sphereDataConfig;
		this->op_sphere = i_op_sphere;
		this->op_sphere_nodealiasing = i_op_sphere_nodealiasing;
		this->setup(i_simVars,
				i_timeSteppersFine, i_timeSteppersCoarse);

		// IMPORTANT: setup initial conditions (inside this->setup) before setting up timesteppers
		// because simulation parameters may change
		this->timeSteppersFine->setup(
					this->simVars->disc.timestepping_method,
					////this->simVars->disc.timestepping_order,
					////this->simVars->disc.timestepping_order2,
					*this->op_sphere[0],
					*this->simVars
				);
		this->timeSteppersCoarse->setup(
					this->simVars_coarse->disc.timestepping_method,
					////this->simVars_coarse->disc.timestepping_order,
					////this->simVars_coarse->disc.timestepping_order2,
					*this->op_sphere[1],
					*this->simVars_coarse
				);

		PInT_Common::setup();
	}
#endif

	void setup(SimulationVariables* i_simVars,
			t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse)
	{

		this->simVars = i_simVars;
		this->simVars_coarse = new SimulationVariables;
		*this->simVars_coarse = *this->simVars;
		this->simVars_coarse->disc.timestepping_method = this->simVars->parareal.coarse_timestepping_method;
		this->simVars_coarse->disc.timestepping_order = this->simVars->parareal.coarse_timestepping_order;
		this->simVars_coarse->disc.timestepping_order2 = this->simVars->parareal.coarse_timestepping_order2;
		this->simVars_coarse->timecontrol.current_timestep_size = this->simVars->parareal.coarse_timestep_size;

		this->timeSteppersFine = i_timeSteppersFine;
		this->timeSteppersCoarse = i_timeSteppersCoarse;


#if SWEET_PARAREAL_SCALAR
		this->timeSteppersFine->setup(*this->simVars);
		this->timeSteppersCoarse->setup(*this->simVars);
#endif

		// All containers contain data defined on the fine spatial grid
		// (including the coarse_data container, with interpolation performed in run_timestep_coarse)
		// Exception: the auxiliary coarse solutions for SL tsm, since they are not summed up with fine solutions
		this->parareal_data_start                      = this->create_new_data_container("fine");
		this->parareal_data_fine                       = this->create_new_data_container("fine");
		this->parareal_data_coarse                     = this->create_new_data_container("fine");
		this->parareal_data_coarse_coarse_mesh         = this->create_new_data_container("coarse");
		this->parareal_data_output                     = this->create_new_data_container("fine");
		this->parareal_data_error                      = this->create_new_data_container("fine");
		this->parareal_data_coarse_previous_timestep   = this->create_new_data_container("coarse");
		this->parareal_data_coarse_previous_time_slice = this->create_new_data_container("coarse");
		this->parareal_data_fine_previous_timestep     = this->create_new_data_container("fine");
		this->parareal_data_fine_previous_time_slice   = this->create_new_data_container("fine");
		this->parareal_data_ref_exact                  = this->create_new_data_container("fine");
		this->parareal_data_fine_exact                 = this->create_new_data_container("fine");
#if SWEET_DEBUG
		this->parareal_data_fine_exact_debug           = this->create_new_data_container("fine");
#endif

		this->parareal_data_debug           = this->create_new_data_container("coarse");

		this->sim_setup_initial_data();
	};


public:

	Parareal_GenericData* create_new_data_container(std::string level)
	{
#if SWEET_PARAREAL_SCALAR
		{
			Parareal_GenericData_Scalar<N>* out = new Parareal_GenericData_Scalar<N>;
			out->allocate_data();
			out->set_time(this->timeframe_end);
			return out;
		}

#elif SWEET_PARAREAL_PLANE
		{
			Parareal_GenericData_PlaneData_Spectral<N>* out = new Parareal_GenericData_PlaneData_Spectral<N>;
			if (level == "fine")
				out->setup_data_config(this->planeDataConfig[0]);
			else if (level == "coarse")
				out->setup_data_config(this->planeDataConfig[1]);
			else
				SWEETError("Wrong level.");
			out->allocate_data();
			return out;
		}

#elif SWEET_PARAREAL_SPHERE
		{
			Parareal_GenericData_SphereData_Spectral<N>* out = new Parareal_GenericData_SphereData_Spectral<N>;
			if (level == "fine")
				out->setup_data_config(this->sphereDataConfig[0]);
			else if (level == "coarse")
				out->setup_data_config(this->sphereDataConfig[1]);
			else
				SWEETError("Wrong level.");
			out->allocate_data();
			return out;
		}
#endif

		// default time set
	}

public:

	/**
	 * Check if the time slice contains an integer number of coarse and fine time steps
	 */
	void sim_check_timesteps(
			double time_slice_size
	)
	{

		// check if seup has been called
		assert(simVars_coarse->timecontrol.current_timestep_size == simVars->parareal.coarse_timestep_size);

		// check if each time slice contains an integer number of fine and coarse time steps
		double eps = 1e-12;
		double mod_coarse = fmod(time_slice_size, simVars_coarse->timecontrol.current_timestep_size);
		double mod_fine = fmod(time_slice_size, simVars->timecontrol.current_timestep_size);
                if ( std::abs(mod_coarse) > eps && std::abs(mod_coarse - time_slice_size) > eps )
			SWEETError("Time slice length must be an integer multiple of the coarse time step! (" + std::to_string(simVars_coarse->timecontrol.current_timestep_size) + ", " + std::to_string(time_slice_size) + ")");
                if ( std::abs(mod_fine) > eps && std::abs(mod_fine - time_slice_size) > eps )
		{
			std::cout << "Number of timesteps: " << this->nb_timesteps_fine << std::endl;
			std::cout << "Mod(dt): " << std::abs(mod_fine) << std::endl;
			std::cout << "Mod(dt) - Dt: " << std::abs(mod_fine - time_slice_size) << std::endl;
			SWEETError("Time slice length must be an integer multiple of the fine time step! (" + std::to_string(simVars->timecontrol.current_timestep_size) + ", " + std::to_string(time_slice_size) + ")");
		}
	};


	/**
	 * Set the start and end of the coarse time step
	 */
	void sim_set_timeframe(
			double i_timeframe_start,	///< start timestamp of coarse time step
			double i_timeframe_end		///< end time stamp of coarse time step
	){
		// check if seup has been called
		assert(simVars_coarse->timecontrol.current_timestep_size == simVars->parareal.coarse_timestep_size);

		if (simVars->parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;
		this->timeframe_start = i_timeframe_start;
		this->timeframe_end = i_timeframe_end;

		this->nb_timesteps_fine = (int)((this->timeframe_end - this->timeframe_start) / simVars->timecontrol.current_timestep_size);
		this->nb_timesteps_coarse = (int)((this->timeframe_end - this->timeframe_start) / simVars_coarse->timecontrol.current_timestep_size);
		if (this->timeframe_start + this->nb_timesteps_fine * simVars->timecontrol.current_timestep_size < this->timeframe_end - 1e-14)
			this->nb_timesteps_fine++;
		if (this->timeframe_start + this->nb_timesteps_coarse * simVars_coarse->timecontrol.current_timestep_size < this->timeframe_end - 1e-14)
			this->nb_timesteps_coarse++;
		assert( std::abs(this->timeframe_start + this->nb_timesteps_fine * simVars->timecontrol.current_timestep_size - this->timeframe_end) < 1e-14);
		assert( std::abs(this->timeframe_start + this->nb_timesteps_coarse * simVars_coarse->timecontrol.current_timestep_size - this->timeframe_end) < 1e-14);

		this->dt_fine = simVars->timecontrol.current_timestep_size;
		this->dt_coarse = simVars->parareal.coarse_timestep_size;

		// set time to parareal_genericdata instances
		this->parareal_data_start->set_time(i_timeframe_end);
		this->parareal_data_fine->set_time(i_timeframe_end);
		this->parareal_data_coarse->set_time(i_timeframe_end);
		this->parareal_data_coarse_coarse_mesh->set_time(i_timeframe_end);
		this->parareal_data_output->set_time(i_timeframe_end);
		this->parareal_data_error->set_time(i_timeframe_end);
		this->parareal_data_coarse_previous_timestep->set_time(i_timeframe_end);
		this->parareal_data_coarse_previous_time_slice->set_time(i_timeframe_end);
		this->parareal_data_fine_previous_timestep->set_time(i_timeframe_end);
		this->parareal_data_fine_previous_time_slice->set_time(i_timeframe_end);
		this->parareal_data_ref_exact->set_time(i_timeframe_end);
		this->parareal_data_fine_exact->set_time(i_timeframe_end);
#if SWEET_DEBUG
		this->parareal_data_fine_exact_debug->set_time(i_timeframe_end);
#endif

	};

	/**
	 * Set the initial data at i_timeframe_start
	 */
	void sim_setup_initial_data()
	{
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_setup_initial_data()" << std::endl;

		///reset();

#if SWEET_PARAREAL_SCALAR
		double u0 = atof(simVars->bogus.var[1].c_str());
		this->parareal_data_start->dataArrays_to_GenericData_Scalar(u0);

#elif SWEET_PARAREAL_PLANE
		PlaneData_Spectral t0_prog_h_pert(planeDataConfig[0]);
		PlaneData_Spectral t0_prog_u(planeDataConfig[0]);
		PlaneData_Spectral t0_prog_v(planeDataConfig[0]);

	#if SWEET_PARAREAL_PLANE_SWE
		SWEPlaneBenchmarksCombined swePlaneBenchmarks;
		swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, *simVars, *op_plane[0]);
	#elif SWEET_PARAREAL_PLANE_BURGERS
		PlaneData_Physical t0_prog_u_phys(t0_prog_u.planeDataConfig[0]);
		PlaneData_Physical t0_prog_v_phys(t0_prog_v.planeDataConfig[0]);
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
	#endif


		this->parareal_data_start->dataArrays_to_GenericData_PlaneData_Spectral(
	#if SWEET_PARAREAL_PLANE_SWE
											t0_prog_h_pert,
	#endif
											t0_prog_u,
											t0_prog_v);


		if (this->simVars->parareal.spatial_coarsening)
			this->parareal_data_coarse_previous_time_slice->restrict(*this->parareal_data_start);
		else
			this->parareal_data_coarse_previous_time_slice->dataArrays_to_GenericData_PlaneData_Spectral(
		#if SWEET_PARAREAL_PLANE_SWE
												t0_prog_h_pert,
		#endif
												t0_prog_u,
												t0_prog_v);

			this->parareal_data_fine_previous_time_slice->dataArrays_to_GenericData_PlaneData_Spectral(
	#if SWEET_PARAREAL_PLANE_SWE
											t0_prog_h_pert,
	#endif
											t0_prog_u,
											t0_prog_v);

#elif SWEET_PARAREAL_SPHERE
		SphereData_Spectral t0_prog_phi_pert(sphereDataConfig[0]);
		SphereData_Spectral t0_prog_vrt(sphereDataConfig[0]);
		SphereData_Spectral t0_prog_div(sphereDataConfig[0]);

		BenchmarksSphereSWE sphereBenchmarks;
		sphereBenchmarks.setup(*simVars, *op_sphere[0]);
		sphereBenchmarks.master->get_initial_state(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);

		this->parareal_data_start->dataArrays_to_GenericData_SphereData_Spectral(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
		if (this->simVars->parareal.spatial_coarsening)
			this->parareal_data_coarse_previous_time_slice->restrict(*this->parareal_data_start);
		else
			this->parareal_data_coarse_previous_time_slice->dataArrays_to_GenericData_SphereData_Spectral(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
		this->parareal_data_fine_previous_time_slice->dataArrays_to_GenericData_SphereData_Spectral(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
#endif
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

		parareal_data_start->set_time(this->timeframe_end);
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

		parareal_data_coarse_previous_time_slice->set_time(this->timeframe_end);
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

		parareal_data_fine_previous_time_slice->set_time(this->timeframe_end);
	};

#if SWEET_PARAREAL == 2
	void sim_set_data_coarse(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_set_data_coarse()" << std::endl;

		// copy to buffers
		*parareal_data_coarse = i_pararealData;

		parareal_data_coarse->set_time(this->timeframe_end);
	};

#endif

#if SWEET_DEBUG
	/**
	* Store exact solution (full fine simulation) at the end of the time slice
	*/
	void sim_set_data_fine_exact(
			Parareal_GenericData &i_pararealData
	)
	{
		*parareal_data_fine_exact_debug = i_pararealData;

		parareal_data_fine_exact_debug->set_time(this->timeframe_end);
	};

	/**
	* Check if solution at time k (end of time slice k-1) is exact (= fine) at iteration k
	*/
	void compare_to_fine_exact(
	)
	{
		Parareal_GenericData* diff = this->create_new_data_container("fine");
		*diff = *parareal_data_fine_exact_debug;
		*diff -= *parareal_data_output;
		///parareal_data_fine_exact_debug->physical_print();
		///parareal_data_output->physical_print();
		std::cout << "Max serial:" << parareal_data_fine_exact_debug->spectral_reduce_maxAbs() << std::endl;
		std::cout << "Max parareal:" << parareal_data_output->spectral_reduce_maxAbs() << std::endl;
		std::cout << "DIFF: " << diff->spectral_reduce_maxAbs() << std::endl;
		assert(diff->spectral_reduce_maxAbs() < 1e-10);
		delete diff;
	};

#endif

	void set_previous_solution(
				std::string tsm_level
				)
	{

		if (tsm_level == "fine")
			timeSteppersFine->master->set_previous_solution(this->parareal_data_fine_previous_time_slice);
		else if (tsm_level == "coarse")
			timeSteppersCoarse->master->set_previous_solution(this->parareal_data_coarse_previous_time_slice);
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
			Parareal_GenericData* io_data,
			std::string tsm_level
	)
	{

		if (tsm_level == "fine")
			timeSteppersFine->master->run_timestep(
						io_data,
						simVars->timecontrol.current_timestep_size,
						simVars->timecontrol.current_simulation_time
					);
		else if (tsm_level == "coarse")
			timeSteppersCoarse->master->run_timestep(
						io_data,
						simVars_coarse->timecontrol.current_timestep_size,
						simVars_coarse->timecontrol.current_simulation_time
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
		simVars->timecontrol.current_timestep_size = this->dt_fine;

		//std::cout << simVars->disc.timestepping_method << std::endl;
		///std::cout << this->SL_tsm.size() << std::endl;
		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if (std::find(this->SL_tsm.begin(), this->SL_tsm.end(), simVars->disc.timestepping_method) != this->SL_tsm.end())
		{
			this->set_previous_solution("fine");
		}

		*(this->parareal_data_fine) = *(this->parareal_data_start);

		int nb_timesteps = 0;
		while (nb_timesteps != this->nb_timesteps_fine)
		{
			// store previous time step
			// to be used as n-1 in SL in the next time slice
			*(this->parareal_data_fine_previous_timestep) = *(this->parareal_data_fine);

			this->run_timestep(this->parareal_data_fine, "fine");

			simVars->timecontrol.current_simulation_time += simVars->timecontrol.current_timestep_size;
			assert(simVars->timecontrol.current_simulation_time <= timeframe_end + 1e-14);
			nb_timesteps++;
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

		simVars_coarse->timecontrol.current_simulation_time = timeframe_start;
		simVars_coarse->timecontrol.max_simulation_time = timeframe_end;
		simVars_coarse->timecontrol.current_timestep_nr = 0;

		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if (std::find(this->SL_tsm.begin(), this->SL_tsm.end(), simVars_coarse->disc.timestepping_method) != this->SL_tsm.end())
			this->set_previous_solution("coarse");

		*(this->parareal_data_coarse) = *(this->parareal_data_start);

		// interpolate to coarse spatial mesh
		if (this->simVars->parareal.spatial_coarsening)
			this->parareal_data_coarse_coarse_mesh->restrict(*this->parareal_data_coarse);
		else
			*this->parareal_data_coarse_coarse_mesh = *this->parareal_data_coarse;

		int nb_timesteps = 0;
		while (nb_timesteps != this->nb_timesteps_coarse)
		{
			// store previous time step
			// to be used as n-1 in SL in the next time slice
			*(this->parareal_data_coarse_previous_timestep) = *(this->parareal_data_coarse_coarse_mesh);

			this->run_timestep(this->parareal_data_coarse, "coarse");
			simVars_coarse->timecontrol.current_simulation_time += simVars_coarse->timecontrol.current_timestep_size;
			assert(simVars_coarse->timecontrol.current_simulation_time <= timeframe_end +  1e-14);
			nb_timesteps++;
		}

		// interpolate to coarse spatial mesh
		if (this->simVars->parareal.spatial_coarsening)
			this->parareal_data_coarse->pad_zeros(*this->parareal_data_coarse_coarse_mesh);
		else
			*this->parareal_data_coarse = *this->parareal_data_coarse_coarse_mesh;

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
	 * return the difference between fine and coarse solution:
	 */
	Parareal_GenericData& get_reference_to_data_timestep_diff()
	{
		return *(this->parareal_data_error);
	};

	void sim_set_data_diff(
			Parareal_GenericData &i_pararealData
	){
		if (simVars->parareal.verbosity > 2)
			std::cout << "sim_set_data_diff()" << std::endl;

		// copy to buffers
		*parareal_data_error = i_pararealData;

		parareal_data_error->set_time(this->timeframe_end);
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
		Parareal_GenericData* tmp = this->create_new_data_container("fine");
		*tmp = *(this->parareal_data_coarse);
                *tmp += *(this->parareal_data_error);

		// compute difference w.r.t. previous output data
		Parareal_GenericData* tmp2 = this->create_new_data_container("fine");
		*tmp2 = *(this->parareal_data_output);
		*tmp2 -= *tmp;
		convergence = tmp2->spectral_reduce_maxAbs();

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


	void output_data_file(
				int iteration_id,
				int time_slice_id,
				bool output_initial_data = false
	)
	{
		if (output_initial_data)
			PInT_Common::output_data_file(
						this->parareal_data_start,
						iteration_id,
						time_slice_id - 1,
						this->timeframe_start
			);
		else
			PInT_Common::output_data_file(
						this->parareal_data_output,
						iteration_id,
						time_slice_id,
						this->timeframe_end
			);

	}

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

		int nvar = N;

		PInT_Common::store_parareal_error(
						this->parareal_data_output,
						parareal_data_ref,
						nvar,
						iteration_id,
						time_slice_id,
						this->timeframe_end,
						path_ref,
						base_solution,
						i_precision
					);

	}



	void output_data_console(
			int iteration_id,
			int time_slice_id
	)
	{
	};


	void check_for_nan_parareal()
	{
		if (this->parareal_data_output->check_for_nan())
			SWEETError("Instability detected in parareal!");
	};

	void delete_data_container(
					Parareal_GenericData* i_data
	)
	{
		if (i_data)
		{
			delete i_data;
			i_data = nullptr;
		}
	}

	~Parareal_SimulationInstance()
	{

			this->delete_data_container(this->parareal_data_start);
			this->delete_data_container(this->parareal_data_fine);
			this->delete_data_container(this->parareal_data_coarse);
			this->delete_data_container(this->parareal_data_coarse_coarse_mesh);
			this->delete_data_container(this->parareal_data_output);
			this->delete_data_container(this->parareal_data_error);
			this->delete_data_container(this->parareal_data_coarse_previous_timestep);
			this->delete_data_container(this->parareal_data_coarse_previous_time_slice);
			this->delete_data_container(this->parareal_data_fine_previous_timestep);
			this->delete_data_container(this->parareal_data_fine_previous_time_slice);
			this->delete_data_container(this->parareal_data_ref_exact);
			this->delete_data_container(this->parareal_data_fine_exact);
	#if SWEET_DEBUG
			this->delete_data_container(this->parareal_data_fine_exact_debug);
	#endif


		if (this->simVars_coarse)
		{
			delete this->simVars_coarse;
			this->simVars_coarse = nullptr;
		}
	}
};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_SIMULATIONINSTANCE_HPP_ */
