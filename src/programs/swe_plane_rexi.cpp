/*
 * SWE with nonlinear part using REXI test
 *
 */

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
#endif

#ifndef SWEET_MPI
	#define  SWEET_MPI 1
#endif


#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>

#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataGridMapping.hpp>
#include <sweet/plane/PlaneDiagnostics.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


#include "swe_plane_rexi/SWE_Plane_TimeSteppers.hpp"



// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;



#ifndef SWEET_MPI
#	define SWEET_MPI 1
#endif

#if SWEET_MPI
#	include <mpi.h>
#endif

#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#if SWEET_PARAREAL
	#include <parareal/Parareal.hpp>
#endif


/// general parameters
SimulationVariables simVars;

double param_initial_freq_x_mul;
double param_initial_freq_y_mul;

class SimulationInstance
#if SWEET_PARAREAL
		:
		public Parareal_SimulationInstance
#endif
{
public:
	// Prognostic variables
	// h: surface height
	// u: velocity in x-direction
	// v: velocity in y-direction
	PlaneData prog_h_pert, prog_u, prog_v;


#if SWEET_GUI
	//visualization variable
	PlaneData vis;
#endif

	// Initial values for comparison with analytical solution
	PlaneData t0_prog_h_pert, t0_prog_u, t0_prog_v;

	// Forcings
	PlaneData force_h_pert, force_u, force_v;

	PlaneDataGridMapping gridMapping;


	// implementation of different time steppers
	SWE_Plane_TimeSteppers timeSteppers;

#if SWEET_PARAREAL
	// implementation of different time steppers
	SWE_Plane_TimeSteppers timeSteppersCoarse;
#endif


	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	class BenchmarkErrors
	{
	public:
		//Max difference to initial conditions
		double t0_error_max_abs_h_pert;
		double t0_error_max_abs_u;
		double t0_error_max_abs_v;

		// Error measures L2 norm
		double analytical_error_rms_h;
		double analytical_error_rms_u;
		double analytical_error_rms_v;

		// Error measures max norm
		double analytical_error_maxabs_h;
		double analytical_error_maxabs_u;
		double analytical_error_maxabs_v;
	};

	BenchmarkErrors benchmark;

	// Finite difference operators
	PlaneOperators op;


	/// Diagnostic measures at initial stage, Initialize with 0
	double diagnostics_energy_start = 0;
	double diagnostics_mass_start = 0;
	double diagnostics_potential_entrophy_start = 0;


	bool compute_error_difference_to_initial_condition = false;
	bool compute_error_to_analytical_solution = false;



public:
	SimulationInstance()	:
		// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_h_pert(planeDataConfig),
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

#if SWEET_GUI
		vis(planeDataConfig),
#endif

		t0_prog_h_pert(planeDataConfig),
		t0_prog_u(planeDataConfig),
		t0_prog_v(planeDataConfig),

		force_h_pert(planeDataConfig),
		force_u(planeDataConfig),
		force_v(planeDataConfig),

		// Initialises operators
		op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
#if SWEET_PARAREAL != 0
		,
		_parareal_data_start_h(planeDataConfig), _parareal_data_start_u(planeDataConfig), _parareal_data_start_v(planeDataConfig),
		_parareal_data_fine_h(planeDataConfig), _parareal_data_fine_u(planeDataConfig), _parareal_data_fine_v(planeDataConfig),
		_parareal_data_coarse_h(planeDataConfig), _parareal_data_coarse_u(planeDataConfig), _parareal_data_coarse_v(planeDataConfig),
		_parareal_data_output_h(planeDataConfig), _parareal_data_output_u(planeDataConfig), _parareal_data_output_v(planeDataConfig),
		_parareal_data_error_h(planeDataConfig), _parareal_data_error_u(planeDataConfig), _parareal_data_error_v(planeDataConfig)
#endif
	{
		// Calls initialisation of the run (e.g. sets u, v, h)
		reset();

#if SWEET_PARAREAL
		parareal_setup();

#endif
	}



	virtual ~SimulationInstance()
	{
	}



	void reset()
	{
		if (simVars.setup.benchmark_scenario_id < 0)
		{
			std::cout << std::endl;
			std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
			SWEPlaneBenchmarks::printScenarioInformation();
			FatalError("Benchmark scenario not selected");
		}


		// Initialise diagnostics
		last_timestep_nr_update_diagnostics = -1;


		benchmark.t0_error_max_abs_h_pert = -1;
		benchmark.t0_error_max_abs_u = -1;
		benchmark.t0_error_max_abs_v = -1;

		benchmark.analytical_error_rms_h = -1;
		benchmark.analytical_error_rms_u = -1;
		benchmark.analytical_error_rms_v = -1;

		benchmark.analytical_error_maxabs_h = -1;
		benchmark.analytical_error_maxabs_u = -1;
		benchmark.analytical_error_maxabs_v = -1;

		simVars.timecontrol.current_timestep_nr = 0;
		simVars.timecontrol.current_simulation_time = 0;

		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		prog_h_pert.spectral_set_all(0, 0);
		prog_u.spectral_set_all(0, 0);
		prog_v.spectral_set_all(0, 0);
#endif

		//Setup prog vars
		prog_h_pert.physical_set_all(simVars.sim.h0);
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);

		//Check if input parameters are adequate for this simulation
		if (simVars.disc.use_staggering && simVars.disc.use_spectral_basis_diffs)
			FatalError("Staggering and spectral basis not supported!");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (simVars.disc.use_staggering ||  !simVars.disc.use_spectral_basis_diffs)
			FatalError("Finite differences and spectral dealisiang should not be used together! Please compile without dealiasing.");
#endif


		if (simVars.disc.use_staggering)
			gridMapping.setup(simVars, planeDataConfig);


		// Waves test case - separate from SWEValidationBench because it depends on certain local input parameters
		auto return_h_perturbed = [] (
				SimulationVariables &i_parameters,
				double x,
				double y
		) -> double
		{
			if (param_initial_freq_x_mul == 0)
				return SWEPlaneBenchmarks::return_h_perturbed(simVars, x, y);

			// Waves scenario
			// Remember to set up initial_freq_x_mul and initial_freq_y_mul
			double dx = x/i_parameters.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
			double dy = y/i_parameters.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
			return std::sin(2.0*dx)*std::cos(2.0*dy) - (1.0/5.0)*std::cos(2.0*dx)*std::sin(4.0*dy);
		};


		auto return_u = [] (
				SimulationVariables &i_parameters,
				double x,
				double y
		) -> double
		{
			if (param_initial_freq_x_mul == 0)
				return SWEPlaneBenchmarks::return_u(simVars, x, y);

			double dx = x/i_parameters.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
			double dy = y/i_parameters.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
			return std::cos(4.0*dx)*std::cos(2.0*dy);
		};


		auto return_v = [] (
				SimulationVariables &i_parameters,
				double x,
				double y
		) -> double
		{
			if (param_initial_freq_x_mul == 0)
				return SWEPlaneBenchmarks::return_v(simVars, x, y);

			double dx = x/i_parameters.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
			double dy = y/i_parameters.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
			return std::cos(2.0*dx)*std::cos(4.0*dy);
		};

		// Set initial conditions given from SWEValidationBenchmarks
		for (int j = 0; j < simVars.disc.res_physical[1]; j++)
		{
			for (int i = 0; i < simVars.disc.res_physical[0]; i++)
			{
				if (simVars.disc.use_staggering) // C-grid
				{
					{
						// h - lives in the center of the cell
						double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

						prog_h_pert.p_physical_set(j, i, return_h_perturbed(simVars, x, y));
						t0_prog_h_pert.p_physical_set(j, i, return_h_perturbed(simVars, x, y));
						force_h_pert.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_h_perturbed(simVars, x, y));
					}


					{
						// u space
						double x = (((double)i)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

						prog_u.p_physical_set(j,i, return_u(simVars, x, y));
						t0_prog_u.p_physical_set(j, i, return_u(simVars, x, y));
						force_u.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_u(simVars, x, y));
					}

					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
						double y = (((double)j)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

						prog_v.p_physical_set(j, i, return_v(simVars, x, y));
						t0_prog_v.p_physical_set(j, i, return_v(simVars, x, y));
						force_v.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_v(simVars, x, y));
					}
				}
				else // A-Grid (colocated grid)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

					prog_h_pert.p_physical_set(j, i, return_h_perturbed(simVars, x, y));
					prog_u.p_physical_set(j, i, return_u(simVars, x, y));
					prog_v.p_physical_set(j, i, return_v(simVars, x, y));

					t0_prog_h_pert.p_physical_set(j, i, return_h_perturbed(simVars, x, y));
					t0_prog_u.p_physical_set(j, i, return_u(simVars, x, y));
					t0_prog_v.p_physical_set(j, i, return_v(simVars, x, y));

					force_h_pert.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_h_perturbed(simVars, x, y));
					force_u.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_u(simVars, x, y));
					force_v.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_v(simVars, x, y));
				}
			}
		}

		// Load data, if requested
		if (simVars.setup.input_data_filenames.size() > 0)
			prog_h_pert.file_physical_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 1)
			prog_u.file_physical_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 2)
			prog_v.file_physical_loadData(simVars.setup.input_data_filenames[2].c_str(), simVars.setup.input_data_binary);

		timeSteppers.setup(
				simVars.disc.timestepping_method,
				simVars.disc.timestepping_order,
				simVars.disc.timestepping_order2,
				op,
				simVars
			);

		if (simVars.misc.compute_errors)
		{
			//Compute difference to initial condition (makes more sense in steady state cases, but useful in others too)
			compute_error_difference_to_initial_condition = true;
					//simVars.setup.benchmark_scenario_id == 2 ||
					//simVars.setup.benchmark_scenario_id == 3 ||
					//simVars.setup.benchmark_scenario_id == 14;

			//Compute difference to analytical solution (makes more sense in linear cases, but might be useful in others too)
			compute_error_to_analytical_solution = timeSteppers.linear_only;
		}
		else
		{
			compute_error_difference_to_initial_condition = false;
			compute_error_to_analytical_solution = false;
		}

		update_diagnostics();

		diagnostics_energy_start = simVars.diag.total_energy;
		diagnostics_mass_start = simVars.diag.total_mass;
		diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;


		timestep_output();
	}


	//Calculate the model diagnostics
	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;


		if (simVars.disc.use_staggering)
		{
			PlaneDiagnostics::update_staggered_huv_to_mass_energy_enstrophy(
					op,
					prog_h_pert,
					prog_u,
					prog_v,
					simVars
			);
		}
		else
		{
			PlaneDiagnostics::update_nonstaggered_huv_to_mass_energy_enstrophy(
					op,
					prog_h_pert,
					prog_u,
					prog_v,
					simVars
			);
		}
	}



	void normal_mode_analysis()
	{
		// dummy time step to get time step size
		if (simVars.sim.CFL >= 0)
			FatalError("Normal mode analysis requires setting fixed time step size");

		simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;

		//run_timestep();

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		const char* filename;
		char buffer_real[1024];

		if (simVars.misc.output_file_name_prefix == "")
			filename = "output_%s_normalmodes.csv";
		else
			filename = simVars.misc.output_file_name_prefix.c_str();


		sprintf(buffer_real, filename, "normal_modes_physical", simVars.timecontrol.current_timestep_size*simVars.misc.output_time_scale);
		std::ofstream file(buffer_real, std::ios_base::trunc);
		std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

		std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

		// use very high precision
		file << std::setprecision(20);

		PlaneData* prog[3] = {&prog_h_pert, &prog_u, &prog_v};

		/*
		 * Maximum number of prognostic variables
		 *
		 * Advection e.g. has only one
		 */
		int max_prog_id = -1;
		if (simVars.pde.id == 0 || simVars.pde.id == 1)
		{
			max_prog_id = 3;
			prog_h_pert.physical_set_zero();
			prog_u.physical_set_zero();
			prog_v.physical_set_zero();
		}
		else
		{
			max_prog_id = 1;
			prog_h_pert.physical_set_zero();
		}

#if 0
		if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
		{
			FatalError("Not yet tested and supported");
			std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
			std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
			std::cout << "Therefore, we also halve the time step size here" << std::endl;

			simVars.timecontrol.current_timestep_size = 0.5*simVars.sim.CFL;
			simVars.sim.CFL = -simVars.timecontrol.current_timestep_size;
		}
#endif

		int num_timesteps = 1;
		if (simVars.disc.normal_mode_analysis_generation >= 10)
		{
			if (simVars.timecontrol.max_timesteps_nr > 0)
				num_timesteps = simVars.timecontrol.max_timesteps_nr;
		}

		if (simVars.timecontrol.max_simulation_time > 0)
			file << "# t " << simVars.timecontrol.max_simulation_time << std::endl;
		else
			file << "# t " << (num_timesteps*(-simVars.sim.CFL)) << std::endl;

		file << "# g " << simVars.sim.gravitation << std::endl;
		file << "# h " << simVars.sim.h0 << std::endl;
		file << "# r " << simVars.sim.earth_radius << std::endl;
		file << "# f " << simVars.sim.coriolis_omega << std::endl;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		int specmodes = planeDataConfig->get_spectral_iteration_range_area(0)+planeDataConfig->get_spectral_iteration_range_area(1);
		file << "# specnummodes " << specmodes << std::endl;
		file << "# specrealresx " << planeDataConfig->spectral_real_modes[0] << std::endl;
		file << "# specrealresy " << planeDataConfig->spectral_real_modes[1] << std::endl;
#endif

		file << "# physresx " << planeDataConfig->physical_res[0] << std::endl;
		file << "# physresy " << planeDataConfig->physical_res[1] << std::endl;
		file << "# normalmodegeneration " << simVars.disc.normal_mode_analysis_generation << std::endl;
		file << "# antialiasing ";

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		file << 1;
#else
		file << 0;
#endif

		file << std::endl;


		// iterate over all prognostic variables
		for (int outer_prog_id = 0; outer_prog_id < max_prog_id; outer_prog_id++)
		{
			if (simVars.disc.normal_mode_analysis_generation == 1 || simVars.disc.normal_mode_analysis_generation == 11)
			{
				// iterate over physical space
				for (std::size_t outer_i = 0; outer_i < planeDataConfig->physical_array_data_number_of_elements; outer_i++)
				{
					// reset time control
					simVars.timecontrol.current_timestep_nr = 0;
					simVars.timecontrol.current_simulation_time = 0;

					std::cout << "." << std::flush;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->physical_set_zero();

					// activate mode
					prog[outer_prog_id]->request_data_physical();
					prog[outer_prog_id]->physical_space_data[outer_i] = 1;

					/*
					 * RUN timestep
					 */

					run_timestep();

					if (simVars.disc.normal_mode_analysis_generation == 1)
					{
						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_physical();
						prog[outer_prog_id]->physical_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							(*prog[inner_prog_id]) /= simVars.timecontrol.current_timestep_size;
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_physical();
						for (std::size_t k = 0; k < planeDataConfig->physical_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->physical_space_data[k];
							if (inner_prog_id != max_prog_id-1 || k != planeDataConfig->physical_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
#if 1
			else if (simVars.disc.normal_mode_analysis_generation == 3 || simVars.disc.normal_mode_analysis_generation == 13)
			{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
				FatalError("Only available with if plane spectral space is activated during compile time!");
#else

				// iterate over spectral space
				for (int r = 0; r < 2; r++)
				{

					for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
					{
						for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
						{
							// reset time control
							simVars.timecontrol.current_timestep_nr = 0;
							simVars.timecontrol.current_simulation_time = 0;

							std::cout << "." << std::flush;

							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
								prog[inner_prog_id]->spectral_set_zero();

							// activate mode via real coefficient
							prog[outer_prog_id]->p_spectral_set(j, i, 1.0);

							/*
							 * RUN timestep
							 */
							run_timestep();


							if (simVars.disc.normal_mode_analysis_generation == 3)
							{
								/*
								 * compute
								 * 1/dt * (U(t+1) - U(t))
								 */
								prog[outer_prog_id]->request_data_spectral();

								std::complex<double> val = prog[outer_prog_id]->p_spectral_get(j, i);
								val = val - 1.0;
								prog[outer_prog_id]->p_spectral_set(j, i, val);

								for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
									(*prog[inner_prog_id]) /= simVars.timecontrol.current_timestep_size;
							}


							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							{
								prog[inner_prog_id]->request_data_spectral();

								/*
								 * REAL
								 */

								for (int r = 0; r < 2; r++)
								{
									for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
									{
										for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
										{
											file << prog[inner_prog_id]->p_spectral_get(j, i).real();
											file << "\t";
										}
									}
								}
							}


							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							{
								/*
								 * IMAG
								 */
								int c = 0;
								for (int r = 0; r < 2; r++)
								{
									for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
									{
										for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
										{
											file << prog[inner_prog_id]->p_spectral_get(j, i).imag();

											if (inner_prog_id != max_prog_id-1 || c != specmodes-1)
												file << "\t";
											else
												file << std::endl;

											c++;
										}
									}
								}
							}
						}
					}
				}
#endif
			}
#else
			else if (simVars.disc.normal_mode_analysis_generation == 3 || simVars.disc.normal_mode_analysis_generation == 13)
			{
				PlaneDataComplex t1(planeDataConfig);
				PlaneDataComplex t2(planeDataConfig);
				PlaneDataComplex t3(planeDataConfig);
				PlaneDataComplex* prog_cplx[3] = {&t1, &t2, &t3};

				// iterate over spectral space
				for (std::size_t outer_i = 0; outer_i < planeDataConfig->spectral_complex_array_data_number_of_elements; outer_i++)
				{
					// reset time control
					simVars.timecontrol.current_timestep_nr = 0;
					simVars.timecontrol.current_simulation_time = 0;

					std::cout << "." << std::flush;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog_cplx[inner_prog_id]->spectral_set_zero();

					// activate mode via real coefficient
					prog_cplx[outer_prog_id]->request_data_spectral();
					prog_cplx[outer_prog_id]->spectral_space_data[outer_i].real(1);

					// convert PlaneDataComplex to PlaneData
					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						*prog[inner_prog_id] = Convert_PlaneDataComplex_To_PlaneData::physical_convert(*prog_cplx[inner_prog_id]);
						prog[inner_prog_id]->spectral_zeroAliasingModes();
					}


					/*
					 * RUN timestep
					 */
					run_timestep();

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->spectral_zeroAliasingModes();
#warning "update this physical_convert maybe to spectral_convert"

						*prog_cplx[inner_prog_id] = Convert_PlaneData_To_PlaneDataComplex::physical_convert(*prog[inner_prog_id]);

						prog_cplx[inner_prog_id]->request_data_spectral();
					}

					if (simVars.disc.normal_mode_analysis_generation == 3)
					{
						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog_cplx[outer_prog_id]->request_data_spectral();
						prog_cplx[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog_cplx[inner_prog_id]->operator*=(1.0/simVars.timecontrol.current_timestep_size);
					}


					// convert PlaneDataComplex to PlaneData
					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog_cplx[inner_prog_id]->request_data_spectral();

						/*
						 * REAL
						 */
						for (std::size_t k = 0; k < planeDataConfig->spectral_complex_array_data_number_of_elements; k++)
						{
							file << prog_cplx[inner_prog_id]->spectral_space_data[k].real();
							file << "\t";
						}

						/*
						 * IMAG
						 */
						for (std::size_t k = 0; k < planeDataConfig->spectral_complex_array_data_number_of_elements; k++)
						{
							file << prog_cplx[inner_prog_id]->spectral_space_data[k].imag();

							if (inner_prog_id != max_prog_id-1 || k != planeDataConfig->spectral_complex_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
#endif
		}
	}


	/**
	 * Execute a single simulation time step
	 */
	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_h_pert, prog_u, prog_v,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);

		// Apply viscosity at posteriori, for all methods explicit diffusion for non spectral schemes and implicit for spectral
		if (simVars.sim.viscosity != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE //TODO: this needs checking

			double dt = simVars.timecontrol.current_timestep_size;

			prog_u = prog_u + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_v = prog_v + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_h_pert = prog_h_pert + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_h_pert, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_h_pert, simVars.sim.viscosity_order)*simVars.sim.viscosity;
#else
			prog_u = op.implicit_diffusion(prog_u, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
			prog_v = op.implicit_diffusion(prog_v, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
			prog_h_pert = op.implicit_diffusion(prog_h_pert, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
#endif
		}

		// advance time step and provide information to parameters
		simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			FatalError("Max simulation time exceeded!");

		timestep_output();
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file(
			const PlaneData &i_planeData,
			const char* i_name	///< name of output variable
		)
	{
		char buffer[1024];

		const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
		i_planeData.file_physical_saveData_ascii(buffer);

		return buffer;
	}



public:
	bool timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (simVars.disc.normal_mode_analysis_generation > 0)
			return false;

		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds-simVars.misc.output_next_sim_seconds*(1e-12) > simVars.timecontrol.current_simulation_time)
			return false;

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		{
			// For output, variables need to be on unstaggered A-grid
			PlaneData t_h(planeDataConfig);
			PlaneData t_u(planeDataConfig);
			PlaneData t_v(planeDataConfig);

			if (simVars.disc.use_staggering) // Remap in case of C-grid
			{
				t_h = prog_h_pert;
				gridMapping.mapCtoA_u(prog_u, t_u);
				gridMapping.mapCtoA_v(prog_v, t_v);
			}
			else
			{
				t_h = prog_h_pert;
				t_u = prog_u;
				t_v = prog_v;
			}


			//Dump  data in csv, if requested
			if (simVars.misc.output_file_name_prefix.size() > 0)
			{
				write_file(t_h, "prog_h_pert");
				write_file(t_u, "prog_u");
				write_file(t_v, "prog_v");

				write_file(op.diff_c_x(prog_v) - op.diff_c_y(prog_u), "prog_q");
			}
		}


		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();
			compute_errors();

			std::stringstream header;
			std::stringstream rows;

			rows << std::setprecision(8);

			// Prefix
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "DATA";
			rows << "DATA";

			// Time
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "\tT";
			rows << "\t" << simVars.timecontrol.current_simulation_time;

#if 1
			// Mass, Energy, Enstrophy
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";
			rows << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;
#endif


#if 1
			if (compute_error_difference_to_initial_condition)
			{
				// Difference to initial condition
				if (simVars.timecontrol.current_timestep_nr == 0)
					header << "\tDIFF_MAXABS_H0\tDIFF_MAXABS_U0\tDIFF_MAXABS_V0";

				rows << "\t" << benchmark.t0_error_max_abs_h_pert << "\t" << benchmark.t0_error_max_abs_u << "\t" << benchmark.t0_error_max_abs_v;
			}
#endif

#if 1
			if (compute_error_to_analytical_solution)
			{
				if (simVars.timecontrol.current_timestep_nr == 0)
					header << "\tREF_DIFF_MAX_H\tREF_DIFF_MAX_U\tREF_DIFF_MAX_V";

				rows << "\t" << benchmark.analytical_error_maxabs_h << "\t" << benchmark.analytical_error_maxabs_u << "\t" << benchmark.analytical_error_maxabs_v;
			}
#endif

			if (simVars.timecontrol.current_timestep_nr == 0)
				o_ostream << header.str() << std::endl;

			o_ostream << rows.str() << std::endl;

#if 1
			if (diagnostics_mass_start > 0.00001 && std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) > 10000000.0)
			{
				std::cerr << "\n DIAGNOSTICS MASS DIFF TOO LARGE:\t" << std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
			}
#endif

		}

		if (simVars.misc.output_each_sim_seconds > 0)
		{
			if (simVars.misc.output_next_sim_seconds == simVars.timecontrol.max_simulation_time)
			{
				simVars.misc.output_next_sim_seconds = std::numeric_limits<double>::infinity();
			}
			else
			{
				while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
					simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;

				if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.max_simulation_time)
					simVars.misc.output_next_sim_seconds = simVars.timecontrol.max_simulation_time;
			}
		}

		return true;
	}


public:
	void compute_errors()
	{
		if (compute_error_difference_to_initial_condition || compute_error_to_analytical_solution)
		{
			/**
			 * Compute difference to initial condition
			 */
			if (compute_error_difference_to_initial_condition)
			{
				benchmark.t0_error_max_abs_h_pert = (prog_h_pert - t0_prog_h_pert).reduce_maxAbs();
				benchmark.t0_error_max_abs_u = (prog_u - t0_prog_u).reduce_maxAbs();
				benchmark.t0_error_max_abs_v = (prog_v - t0_prog_v).reduce_maxAbs();
			}

			// Calculate linear exact solution, if compute error requests
			if (compute_error_to_analytical_solution)
			{
				// Analytical solution at specific time on A-grid
				PlaneData ts_h_pert = t0_prog_h_pert;
				PlaneData ts_u = t0_prog_u;
				PlaneData ts_v = t0_prog_v;

				// Run exact solution for linear case
				timeSteppers.l_direct->run_timestep(
						ts_h_pert, ts_u, ts_v,
						simVars.timecontrol.current_simulation_time,
						0			// initial condition given at time 0
				);

				benchmark.analytical_error_rms_h = (ts_h_pert-prog_h_pert).reduce_rms_quad();
				benchmark.analytical_error_rms_u = (ts_u-prog_u).reduce_rms_quad();
				benchmark.analytical_error_rms_v = (ts_v-prog_v).reduce_rms_quad();

				benchmark.analytical_error_maxabs_h = (ts_h_pert-prog_h_pert).reduce_maxAbs();
				benchmark.analytical_error_maxabs_u = (ts_u-prog_u).reduce_maxAbs();
				benchmark.analytical_error_maxabs_v = (ts_v-prog_v).reduce_maxAbs();
			}
		}
	}



public:
	bool should_quit()
	{
		if (
				simVars.timecontrol.max_timesteps_nr != -1 &&
				simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr
		)
			return true;

		if (!std::isinf(simVars.timecontrol.max_simulation_time))
			if (simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time+simVars.timecontrol.max_simulation_time*1e-10)	// care about roundoff errors with 1e-10
				return true;

		return false;
	}


#if SWEET_GUI

	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}


	struct VisStuff
	{
		const PlaneData* data;
		const char *description;
	};

	/**
	 * Arrays for online visualisation and their textual description
	 */
	VisStuff vis_arrays[3] =
	{
			{&prog_h_pert,	"h"},
			{&prog_u,	"u"},
			{&prog_v,	"v"},
	};



	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data
	)
	{
		if (simVars.misc.vis_id < 0)
		{
			// Analytical solution at specific time on A-grid
			PlaneData ts_h_pert = t0_prog_h_pert;
			PlaneData ts_u = t0_prog_u;
			PlaneData ts_v = t0_prog_v;

			// Run exact solution for linear case
			timeSteppers.l_direct->run_timestep(
					ts_h_pert, ts_u, ts_v,
					simVars.timecontrol.current_simulation_time,
					0			// initial condition given at time 0
			);

#if 0
			switch(simVars.misc.vis_id)
			{
			case -1:
				vis = ts_u+simVars.sim.h0;			//Exact solution
				break;

			case -2:
				vis = ts_u-prog_u;	// difference to exact solution
				break;

			case -3:
				vis = t0_prog_u-prog_u;	// difference to initial condition
				break;
			}
#else
			switch(simVars.misc.vis_id)
			{
			case -1:
				vis = ts_h_pert+simVars.sim.h0;			//Exact solution
//				vis = ts_u;
				break;

			case -2:
				vis = ts_h_pert-prog_h_pert;	// difference to exact solution
				break;

			case -3:
				vis = t0_prog_h_pert-prog_h_pert;	// difference to initial condition
				break;
			}
#endif
			*o_dataArray = &vis;
			*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
			return;
		}

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		if (id == 0)
		{
			vis = *vis_arrays[id].data+simVars.sim.h0;
			*o_dataArray = &vis;
		}
		else
			*o_dataArray = vis_arrays[id].data;

		vis=**o_dataArray;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostic values if required
		update_diagnostics();

		const char* description = "";
		if (simVars.misc.vis_id >= 0)
		{
			int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
			description = vis_arrays[id].description;
		}
		else
		{
			switch (simVars.misc.vis_id)
			{
			case -1:
				description = "Direct solution for h (linear only)";
				break;


			case -2:
				description = "Diff in h to exact linear solution";
				break;

			case -3:
				description = "Diff in h to initial condition";
				break;
			}
		}

		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string, "Time (days): %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy,
				vis.reduce_max(),
				vis.reduce_min() );

		return title_string;
	}



	void vis_pause()
	{
		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			simVars.misc.vis_id++;
			break;

		case 'V':
			simVars.misc.vis_id--;
			break;

		case 'c':
			// dump data arrays
			prog_h_pert.file_physical_saveData_ascii("swe_rexi_dump_h.csv");
			prog_u.file_physical_saveData_ascii("swe_rexi_dump_u.csv");
			prog_v.file_physical_saveData_ascii("swe_rexi_dump_v.csv");
			break;

		case 'C':
			// dump data arrays to VTK
			prog_h_pert.file_physical_saveData_vtk("swe_rexi_dump_h.vtk", "Height");
			prog_u.file_physical_saveData_vtk("swe_rexi_dump_u.vtk", "U-Velocity");
			prog_v.file_physical_saveData_vtk("swe_rexi_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_h_pert.file_physical_loadData("swe_rexi_dump_h.csv", simVars.setup.input_data_binary);
			prog_u.file_physical_loadData("swe_rexi_dump_u.csv", simVars.setup.input_data_binary);
			prog_v.file_physical_loadData("swe_rexi_dump_v.csv", simVars.setup.input_data_binary);
			break;
		}
	}
#endif


	bool instability_detected()
	{
		return !(	prog_h_pert.reduce_boolean_all_finite() &&
					prog_u.reduce_boolean_all_finite() &&
					prog_v.reduce_boolean_all_finite()
				);
	}


#if SWEET_PARAREAL

	/******************************************************
	 ******************************************************
	 *       ************** PARAREAL **************
	 ******************************************************
	 ******************************************************/

	PlaneData _parareal_data_start_h, _parareal_data_start_u, _parareal_data_start_v;
	Parareal_Data_PlaneData<3> parareal_data_start;

	PlaneData _parareal_data_fine_h, _parareal_data_fine_u, _parareal_data_fine_v;
	Parareal_Data_PlaneData<3> parareal_data_fine;

	PlaneData _parareal_data_coarse_h, _parareal_data_coarse_u, _parareal_data_coarse_v;
	Parareal_Data_PlaneData<3> parareal_data_coarse;

	PlaneData _parareal_data_output_h, _parareal_data_output_u, _parareal_data_output_v;
	Parareal_Data_PlaneData<3> parareal_data_output;

	PlaneData _parareal_data_error_h, _parareal_data_error_u, _parareal_data_error_v;
	Parareal_Data_PlaneData<3> parareal_data_error;

	double timeframe_start = -1;
	double timeframe_end = -1;

	bool output_data_valid = false;

	void parareal_setup()
	{
		{
			PlaneData* data_array[3] = {&_parareal_data_start_h, &_parareal_data_start_u, &_parareal_data_start_v};
			parareal_data_start.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_fine_h, &_parareal_data_fine_u, &_parareal_data_fine_v};
			parareal_data_fine.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_coarse_h, &_parareal_data_coarse_u, &_parareal_data_coarse_v};
			parareal_data_coarse.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_output_h, &_parareal_data_output_u, &_parareal_data_output_v};
			parareal_data_output.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_error_h, &_parareal_data_error_u, &_parareal_data_error_v};
			parareal_data_error.setup(data_array);
		}

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

		output_data_valid = false;
	}



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

		reset();


		*parareal_data_start.data_arrays[0] = prog_h_pert;
		*parareal_data_start.data_arrays[1] = prog_u;
		*parareal_data_start.data_arrays[2] = prog_v;

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

		// cast to pararealPlaneData stuff
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

		while (simVars.timecontrol.current_simulation_time != timeframe_end)
		{
			run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// copy to buffers
		*parareal_data_fine.data_arrays[0] = prog_h_pert;
		*parareal_data_fine.data_arrays[1] = prog_u;
		*parareal_data_fine.data_arrays[2] = prog_v;
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

		prog_h_pert = *parareal_data_start.data_arrays[0];
		prog_u = *parareal_data_start.data_arrays[1];
		prog_v = *parareal_data_start.data_arrays[2];

		timestepping_swe_plane_rexi.run_timestep_implicit_ts(
				prog_h_pert, prog_u, prog_v,
				timeframe_end - timeframe_start,
				op,
				simVars
		);


		// copy to buffers
		*parareal_data_coarse.data_arrays[0] = prog_h_pert;
		*parareal_data_coarse.data_arrays[1] = prog_u;
		*parareal_data_coarse.data_arrays[2] = prog_v;
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

		for (int k = 0; k < 3; k++)
			*parareal_data_error.data_arrays[k] = *parareal_data_fine.data_arrays[k] - *parareal_data_coarse.data_arrays[k];
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
		double convergence = -1;

		if (!i_compute_convergence_test || !output_data_valid)
		{
			for (int k = 0; k < 3; k++)
				*parareal_data_output.data_arrays[k] = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			output_data_valid = true;
			return convergence;
		}



		for (int k = 0; k < 3; k++)
		{
			PlaneData tmp = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			convergence = std::max(
					convergence,
					(*parareal_data_output.data_arrays[k]-tmp).reduce_maxAbs()
				);

			*parareal_data_output.data_arrays[k] = tmp;
		}

		simVars.timecontrol.current_simulation_time = timeframe_end;
		prog_h_pert = *parareal_data_output.data_arrays[0];
		prog_u = *parareal_data_output.data_arrays[1];
		prog_v = *parareal_data_output.data_arrays[2];

		if (compute_error_to_analytical_solution)
		{
			if (simVars.misc.compute_errors > 0)// && simVars.pde.use_nonlinear_equations == 0)
			{
				compute_errors();
				std::cout << "maxabs error compared to analytical solution: " << benchmark.analytical_error_maxabs_h << std::endl;
			}
		}

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


	void output_data_file(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
		Parareal_Data_PlaneData<3>& data = (Parareal_Data_PlaneData<3>&)i_data;

		std::ostringstream ss;
		ss << "output_iter" << iteration_id << "_slice" << time_slice_id << ".vtk";

		std::string filename = ss.str();

		data.data_arrays[0]->file_physical_saveData_vtk(filename.c_str(), filename.c_str());
	}


	void output_data_console(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
	}

#endif
};



int main(int i_argc, char *i_argv[])
{
#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

#if SWEET_MPI

	#if SWEET_THREADING
		int provided;
		MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE)
			FatalError("MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye...");
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

	MemBlockAlloc::setup();

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"initial-freq-x-mul",		/// frequency multipliers for special scenario setup
			"initial-freq-y-mul",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 0;  //frequency in x for waves test case
	simVars.bogus.var[1] = 0;  //frequency in y for waves test case

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--initial-freq-x-mul [float]" << std::endl;
		std::cout << "	--initial-freq-y-mul [float]" << std::endl;
		std::cout << "" << std::endl;


#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}
	// Frequency for certain initial conditions
	param_initial_freq_x_mul = simVars.bogus.var[0];
	param_initial_freq_y_mul = simVars.bogus.var[1];

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);
	planeDataConfig->printInformation();

	// Print header
	std::cout << std::endl;
	simVars.outputConfig();
	std::cout << "Computing error: " << simVars.misc.compute_errors << std::endl;
	std::cout << std::endl;

	std::ostringstream buf;
	buf << std::setprecision(14);


#if SWEET_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	// only start simulation and time stepping for first rank
	if (mpi_rank == 0)
#endif
	{
#if SWEET_MPI
		std::cout << "MPI_RANK: " << mpi_rank << std::endl;
#endif

#if SWEET_PARAREAL
		if (simVars.parareal.enabled)
		{
			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller_Serial<SimulationInstance> parareal_Controller_Serial;

			// setup controller. This initializes several simulation instances
			parareal_Controller_Serial.setup(&simVars.parareal);

			// execute the simulation
			parareal_Controller_Serial.run();
		}
		else
#endif


#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			VisSweet<SimulationInstance> visSweet(simulationSWE);
			delete simulationSWE;
		}
		else
#endif
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			//Setting initial conditions and workspace - in case there is no GUI

			// also initializes diagnostics
			simulationSWE->reset();

			//Time counter
			Stopwatch time;

#if SWEET_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif

			// Start counting time
			time.reset();


			if (simVars.disc.normal_mode_analysis_generation > 0)
			{
				simulationSWE->normal_mode_analysis();
			}
			else
			{
				// Main time loop
				while (true)
				{
					// Stop simulation if requested
					if (simulationSWE->should_quit())
						break;

					// Main call for timestep run
					simulationSWE->run_timestep();

					// Instability
					if (simulationSWE->instability_detected())
					{
						std::cout << "INSTABILITY DETECTED" << std::endl;
						break;
					}
				}
			}

			// Stop counting time
			time.stop();

			double seconds = time();

			// End of run output results
			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;

			simulationSWE->compute_errors();

			if (simVars.misc.verbosity > 0)
			{
				std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((simVars.diag.total_energy-simulationSWE->diagnostics_energy_start)/simulationSWE->diagnostics_energy_start) << std::endl;
				std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((simVars.diag.total_mass-simulationSWE->diagnostics_mass_start)/simulationSWE->diagnostics_mass_start) << std::endl;
				std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((simVars.diag.total_potential_enstrophy-simulationSWE->diagnostics_potential_entrophy_start)/simulationSWE->diagnostics_potential_entrophy_start) << std::endl;

				if (simVars.misc.compute_errors)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark.t0_error_max_abs_h_pert << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark.t0_error_max_abs_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark.t0_error_max_abs_v << std::endl;
				}
			}


			if (simulationSWE->compute_error_to_analytical_solution)
			{
				std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << simulationSWE->benchmark.analytical_error_rms_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << simulationSWE->benchmark.analytical_error_rms_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << simulationSWE->benchmark.analytical_error_rms_v << std::endl;

				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << simulationSWE->benchmark.analytical_error_maxabs_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << simulationSWE->benchmark.analytical_error_maxabs_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << simulationSWE->benchmark.analytical_error_maxabs_v << std::endl;
			}

			delete simulationSWE;
		}
	}
#if SWEET_MPI
	else
	{
		if (simVars.disc.timestepping_method.find("rexi") != std::string::npos)
		{
			PlaneOperators op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

			SWE_Plane_TS_l_rexi rexiSWE(simVars, op);

			/*
			 * Setup our little dog REXI
			 */
			rexiSWE.setup(simVars.rexi);

			bool run = true;

			PlaneData prog_h_pert(planeDataConfig);
			PlaneData prog_u(planeDataConfig);
			PlaneData prog_v(planeDataConfig);

			MPI_Barrier(MPI_COMM_WORLD);

			do
			{
				// REXI time stepping
				rexiSWE.run_timestep(
						prog_h_pert,
						prog_u,
						prog_v,
						-simVars.sim.CFL,		///< if this value is not equal to 0, use this time step size instead of computing one
						simVars.timecontrol.current_simulation_time

				);
			}
			while(!rexiSWE.final_timestep);
		}
	}

	if (simVars.disc.timestepping_method.find("rexi") != std::string::npos)
	{
		// synchronize REXI
		if (mpi_rank == 0)
			SWE_Plane_TS_l_rexi::MPI_quitWorkers(planeDataConfig);
	}

	MPI_Finalize();
#endif

	return 0;
}

