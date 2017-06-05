/*
 * SWE with nonlinear part using REXI test
 *
 */

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
#endif



#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>

#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDiagnostics.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>
#include <sweet/plane/Staggering.hpp>
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

/// specific parameters

bool param_compute_error;

double param_initial_freq_x_mul;
double param_initial_freq_y_mul;

//Diagnostic measures at initial stage
double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;



class SimulationInstance
#if SWEET_PARAREAL
		:
		public Parareal_SimulationInstance
#endif
{
public:
	// Prognostic variables
	PlaneData prog_h, prog_u, prog_v;

	ScalarDataArray pos_x, pos_y;
	ScalarDataArray posa_x, posa_y;

#if SWEET_GUI
	//visualization variable
	PlaneData vis;
#endif

	// Initial values for comparison with analytical solution
	PlaneData t0_prog_h, t0_prog_u, t0_prog_v;

	// Forcings
	PlaneData force_h, force_u, force_v;


	// Staggering
	Staggering staggering;


	// implementation of different time steppers
	SWE_Plane_TimeSteppers timeSteppers;


	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	class BenchmarkErrors
	{
	public:
		//Max difference to initial conditions
		double t0_error_max_abs_h;
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

	// Interpolation stuff
	PlaneDataSampler sampler2D;

public:
	SimulationInstance()	:
	// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_h(planeDataConfig),
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		pos_x(planeDataConfig->physical_array_data_number_of_elements),
		pos_y(planeDataConfig->physical_array_data_number_of_elements),

		posa_x(planeDataConfig->physical_array_data_number_of_elements),
		posa_y(planeDataConfig->physical_array_data_number_of_elements),

#if SWEET_GUI
		vis(planeDataConfig),
#endif

		t0_prog_h(planeDataConfig),
		t0_prog_u(planeDataConfig),
		t0_prog_v(planeDataConfig),

		force_h(planeDataConfig),
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


		PlaneData tmp_x(op.planeDataConfig);
		tmp_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = ((double)i)*simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
			}
		);
		PlaneData tmp_y(op.planeDataConfig);
		tmp_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = ((double)j)*simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];
			}
		);

		// Initialize arrival points with h position
		pos_x = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_x);
		pos_y = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_y);

		// Initialize arrival points with h position
		posa_x = pos_x+0.5*simVars.disc.cell_size[0];
		posa_y = pos_y+0.5*simVars.disc.cell_size[1];

		// Initialise diagnostics
		last_timestep_nr_update_diagnostics = -1;

		benchmark.t0_error_max_abs_h = 0;
		benchmark.t0_error_max_abs_u = 0;
		benchmark.t0_error_max_abs_v = 0;

		benchmark.analytical_error_rms_h = 0;
		benchmark.analytical_error_rms_u = 0;
		benchmark.analytical_error_rms_v = 0;

		benchmark.analytical_error_maxabs_h = 0;
		benchmark.analytical_error_maxabs_u = 0;
		benchmark.analytical_error_maxabs_v = 0;

		simVars.timecontrol.current_timestep_nr = 0;
		simVars.timecontrol.current_simulation_time = 0;

		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		prog_h.spectral_set_all(0, 0);
		prog_u.spectral_set_all(0, 0);
		prog_v.spectral_set_all(0, 0);
#endif

		//Setup prog vars
		prog_h.physical_set_all(simVars.sim.h0);
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);

		//Check if input parameters are adequate for this simulation
		if (simVars.disc.use_staggering && simVars.disc.use_spectral_basis_diffs)
			FatalError("Staggering and spectral basis not supported!");

		if (simVars.disc.use_staggering &&
				( simVars.disc.timestepping_method != SimulationVariables::Discretization::RUNGE_KUTTA_EXPLICIT &&
					simVars.disc.timestepping_method != SimulationVariables::Discretization::SEMI_LAGRANGIAN_ADVECTION_ONLY
				)
		)
			FatalError("Staggering only supported for standard time stepping mode 0 and 4!");

		if (simVars.disc.use_staggering && param_compute_error)
			std::cerr << "Warning: Staggered data will be interpolated to/from A-grid for exact linear solution" << std::endl;

		if (simVars.pde.use_nonlinear_equations > 0 && param_compute_error)
			FatalError("Exact solution not possible in general for nonlinear SWE");

		if (simVars.disc.use_staggering)
			staggering.setup_c_staggering();
		else
			staggering.setup_a_staggering();

		// Setup sampler for future interpolations
		sampler2D.setup(simVars.sim.domain_size, planeDataConfig);

#if 0		// Setup semi-lag
		semiLagrangian.setup(simVars.sim.domain_size, planeDataConfig);


		PlaneData tmp_x(planeDataConfig);
		tmp_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = ((double)i)*simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
			}
		);
		PlaneData tmp_y(planeDataConfig);
		tmp_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = ((double)j)*simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];
			}
		);

		// Initialize arrival points with h position
		pos_x = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_x);
		pos_y = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_y);

		// Initialize arrival points with h position
		posx_a = pos_x+0.5*simVars.disc.cell_size[0];
		posy_a = pos_y+0.5*simVars.disc.cell_size[1];
#endif

		// Waves test case - separate from SWEValidationBench because it depends on certain local input parameters
		auto return_h = [] (
				SimulationVariables &i_parameters,
				double x,
				double y
		) -> double
		{
			if (param_initial_freq_x_mul == 0)
				return SWEPlaneBenchmarks::return_h(simVars, x, y);

			// Waves scenario
			// Remember to set up initial_freq_x_mul and initial_freq_y_mul
			double dx = x/i_parameters.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
			double dy = y/i_parameters.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
			return std::sin(2.0*dx)*std::cos(2.0*dy) - (1.0/5.0)*std::cos(2.0*dx)*std::sin(4.0*dy) + i_parameters.sim.h0;
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

						prog_h.p_physical_set(j, i, return_h(simVars, x, y));
						t0_prog_h.p_physical_set(j, i, return_h(simVars, x, y));
						force_h.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_h(simVars, x, y));
					}


					{
#if 0
						//Coriolis term - lives in the corner of the cells
						if (simVars.pde.use_nonlinear_equations)
						{
							//PXT: had some -0.5 on i and j (why??)
							double x = (((double)i)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
							double y = (((double)j)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
						}
#endif
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

					prog_h.p_physical_set(j, i, return_h(simVars, x, y));
					prog_u.p_physical_set(j, i, return_u(simVars, x, y));
					prog_v.p_physical_set(j, i, return_v(simVars, x, y));

					t0_prog_h.p_physical_set(j, i, return_h(simVars, x, y));
					t0_prog_u.p_physical_set(j, i, return_u(simVars, x, y));
					t0_prog_v.p_physical_set(j, i, return_v(simVars, x, y));

					force_h.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_h(simVars, x, y));
					force_u.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_u(simVars, x, y));
					force_v.p_physical_set(j, i, SWEPlaneBenchmarks::return_force_v(simVars, x, y));
				}
			}
		}

		// Load data, if requested
		if (simVars.setup.input_data_filenames.size() > 0)
			prog_h.file_physical_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 1)
		{
			prog_u.file_physical_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);
			if (simVars.disc.use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

		if (simVars.setup.input_data_filenames.size() > 2)
		{
			prog_v.file_physical_loadData(simVars.setup.input_data_filenames[2].c_str(), simVars.setup.input_data_binary);
			if (simVars.disc.use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

		timeSteppers.setup(simVars.disc.timestepping_method_string, op, simVars);

		// Print output info (if gui is disabled, this is done in main
		if (simVars.misc.gui_enabled)
			timestep_output();
	}


	//Calculate the model diagnostics
	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;


		PlaneDiagnostics::update_nonstaggered_huv_to_mass_energy_enstrophy(
				op,
				prog_h,
				prog_u,
				prog_v,
				simVars
		);
	}



	void normal_mode_analysis()
	{
		// dummy time step to get time step size
		run_timestep();

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		char buffer_real[1024];
		const char* filename = simVars.misc.output_file_name_prefix.c_str();

		sprintf(buffer_real, filename, "normal_modes_physical", simVars.timecontrol.current_timestep_size*simVars.misc.output_time_scale);
		std::ofstream file(buffer_real, std::ios_base::trunc);
		std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

		std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

		// use very high precision
		file << std::setprecision(20);

		PlaneData* prog[3] = {&prog_h, &prog_u, &prog_v};

		int max_prog_id = -1;

		if (simVars.pde.id == 0 || simVars.pde.id == 1)
		{
			max_prog_id = 3;
			prog_h.physical_set_zero();
			prog_u.physical_set_zero();
			prog_v.physical_set_zero();
		}
		else
		{
			max_prog_id = 1;
			prog_h.physical_set_zero();
		}


		if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
		{
			FatalError("Not yet tested and supported");
			std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
			std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
			std::cout << "Therefore, we also halve the time step size here" << std::endl;

			simVars.timecontrol.current_timestep_size = 0.5*simVars.sim.CFL;
			simVars.sim.CFL = -simVars.timecontrol.current_timestep_size;
		}

		// iterate over all prognostic variables
		for (int outer_prog_id = 0; outer_prog_id < max_prog_id; outer_prog_id++)
		{
			if (simVars.disc.normal_mode_analysis_generation == 1)
			{
				// iterate over physical space
				for (std::size_t outer_i = 0; outer_i < planeDataConfig->physical_array_data_number_of_elements; outer_i++)
				{
					std::cout << "normal mode analysis for prog " << outer_prog_id << ", idx " << outer_i << std::endl;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->physical_set_zero();

					// activate mode
					prog[outer_prog_id]->request_data_physical();
					prog[outer_prog_id]->physical_space_data[outer_i] = 1;

					// In case of a multi-step scheme, reset it!
					if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
					{
						FatalError("Not yet implemented");
						//timestepping_explicit_leapfrog.resetAndSetup(prog_h, simVars.disc.timestepping_order, simVars.disc.leapfrog_robert_asselin_filter);
						run_timestep();
					}

					/*
					 * RUN timestep
					 */

					run_timestep();

					/*
					 * compute
					 * 1/dt * (U(t+1) - U(t))
					 */
					prog[outer_prog_id]->request_data_physical();
					prog[outer_prog_id]->physical_space_data[outer_i] -= 1.0;


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
			else if (simVars.disc.normal_mode_analysis_generation == 2)
			{
				// iterate over physical space
				for (std::size_t outer_i = 0; outer_i < planeDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					for (int imag_i = 0; imag_i < 2; imag_i++)
					{
						std::cout << "normal mode analysis for prog " << outer_prog_id << ", idx " << outer_i << std::endl;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->spectral_set_zero();

						// activate mode
						if (imag_i)
							prog[outer_prog_id]->spectral_space_data[outer_i].imag(1);
						else
							prog[outer_prog_id]->spectral_space_data[outer_i].real(1);

						// In case of a multi-step scheme, reset it!
						if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
						{
							FatalError("Not yet implemented");
							//timestepping_explicit_leapfrog.resetAndSetup(prog_h, simVars.disc.timestepping_order, simVars.disc.leapfrog_robert_asselin_filter);
							run_timestep();
						}

						/*
						 * RUN timestep
						 */
						run_timestep();


						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_spectral();
						prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->operator*=(1.0/simVars.timecontrol.current_timestep_size);

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							prog[inner_prog_id]->request_data_spectral();
							for (std::size_t k = 0; k < planeDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].real();
								file << "\t";
							}
						}

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							prog[inner_prog_id]->request_data_spectral();
							for (std::size_t k = 0; k < planeDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].imag();
								if (inner_prog_id != max_prog_id-1 || k != planeDataConfig->spectral_array_data_number_of_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
			}
			else if (simVars.disc.normal_mode_analysis_generation == 3)
			{
				// iterate over physical space
				for (std::size_t outer_i = 0; outer_i < planeDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					std::cout << "normal mode analysis for prog " << outer_prog_id << ", idx " << outer_i << std::endl;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->spectral_set_zero();

					// activate mode via real coefficient
					prog[outer_prog_id]->request_data_spectral();
					prog[outer_prog_id]->spectral_space_data[outer_i].real(1);

					// In case of a multi-step scheme, reset it!
					if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
					{
						FatalError("Not yet implemented");
						//timestepping_explicit_leapfrog.resetAndSetup(prog_h, simVars.disc.timestepping_order, simVars.disc.leapfrog_robert_asselin_filter);
						run_timestep();
					}

					/*
					 * RUN timestep
					 */
					run_timestep();

					/*
					 * compute
					 * 1/dt * (U(t+1) - U(t))
					 */
					prog[outer_prog_id]->request_data_spectral();
					prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->operator*=(1.0/simVars.timecontrol.current_timestep_size);


					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_spectral();

						for (std::size_t k = 0; k < planeDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].real();
							file << "\t";
						}
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_spectral();

						for (std::size_t k = 0; k < planeDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].imag();

							if (inner_prog_id != max_prog_id-1 || k != planeDataConfig->spectral_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
		}
	}


	/**
	 * Execute a single simulation time step
	 */
	void run_timestep()
	{
		double o_dt;

		simVars.timecontrol.current_timestep_size = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 0);

		timeSteppers.master->run_timestep(
				prog_h, prog_u, prog_v,
				o_dt,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.max_simulation_time
			);

#if 0

		else if (simVars.disc.timestepping_method == SimulationVariables::Discretization::SEMI_LAGRANGIAN_ADVECTION_ONLY)
		{
			// Semi-Lagrangian advection - velocities are kept constant, and
			//  h is transported according to given initial wind
			// *** Valid only for non-divergent velocity fields **//
			// See test case scenario 17 for this ts mode

			assert(simVars.sim.CFL < 0);
			o_dt = -simVars.sim.CFL;

			semiLagrangian.semi_lag_departure_points_settls(
							prog_u_prev,
							prog_v_prev,
							prog_u,
							prog_v,
							posx_a,
							posy_a,
							o_dt,
							posx_d,
							posy_d,
							staggering
					);

			//prog_u_prev = prog_u;
			//prog_v_prev = prog_v;
			prog_h_prev = prog_h;

			//Departure points are set for physical space
			//    interpolate to h grid
			sampler2D.bicubic_scalar(
					prog_h_prev,
					posx_d,
					posy_d,
					prog_h,
					staggering.h[0],
					staggering.h[1]
			);

		}
		else if (simVars.disc.timestepping_method == SimulationVariables::Discretization::SEMI_LAGRANGIAN_SEMI_IMPLICIT)
		{
			// Semi-Lagrangian Semi-implicit Spectral
			assert(simVars.sim.CFL < 0);
			o_dt = -simVars.sim.CFL;

			timestepping_swe_plane_rexi.run_timestep_cn_sl_ts(
					prog_h, prog_u, prog_v,
					prog_h_prev, prog_u_prev, prog_v_prev,
					posx_a,	posy_a,
					o_dt,
					simVars.pde.use_nonlinear_equations,
					simVars,
					op,
					sampler2D,
					semiLagrangian
			);
		}
#endif

		//Apply viscosity at posteriori, for all methods explicit diffusion for non spectral schemes and implicit for spectral
		if (simVars.sim.viscosity != 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE //TODO: this needs checking
			prog_u = prog_u + pow(-1,simVars.sim.viscosity_order/2)* o_dt*op.diffN_x(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*o_dt*op.diffN_y(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_v = prog_v + pow(-1,simVars.sim.viscosity_order/2)* o_dt*op.diffN_x(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*o_dt*op.diffN_y(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_h = prog_h + pow(-1,simVars.sim.viscosity_order/2)* o_dt*op.diffN_x(prog_h, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*o_dt*op.diffN_y(prog_h, simVars.sim.viscosity_order)*simVars.sim.viscosity;
#else
			prog_u = op.implicit_diffusion(prog_u, o_dt*simVars.sim.viscosity, simVars.sim.viscosity_order );
			prog_v = op.implicit_diffusion(prog_v, o_dt*simVars.sim.viscosity, simVars.sim.viscosity_order );
			prog_h = op.implicit_diffusion(prog_h, o_dt*simVars.sim.viscosity, simVars.sim.viscosity_order );
#endif
		}

		// advance time step and provide information to parameters
		simVars.timecontrol.current_timestep_size = o_dt;
		simVars.timecontrol.current_simulation_time += o_dt;
		simVars.timecontrol.current_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			FatalError("Max simulation time exceeded!");
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
		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
			return false;

		//Dump  data in csv, if requested
		if (simVars.misc.output_file_name_prefix.size() > 0)
		{
			write_file(prog_h, "prog_h");
			write_file(prog_u, "prog_u");
			write_file(prog_v, "prog_v");

			write_file(op.diff_c_x(prog_v) - op.diff_c_y(prog_u), "prog_q");
		}


		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();

			// Print header
			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				//if ((simVars.setup.scenario >= 0 && simVars.setup.scenario <= 4) || simVars.setup.scenario == 13)
				o_ostream << "\tDIFF_H0\tDIFF_U0\tDIFF_V0";

				if (param_compute_error && simVars.pde.use_nonlinear_equations==0){
					o_ostream << "\tANAL_DIFF_RMS_P\tANAL_DIFF_RMS_U\tANAL_DIFF_RMS_V";
					o_ostream << "\tANAL_DIFF_MAX_P\tANAL_DIFF_MAX_U\tANAL_DIFF_MAX_V";
				}
				o_ostream << std::endl;
			}

			//Print simulation time, energy and pot enstrophy
			o_ostream << std::setprecision(8) << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;


			// PXT: I didn't know where to put this to work with and without GUI - if removed crashes when gui=enable
			if (diagnostics_mass_start==0)
				diagnostics_mass_start=simVars.diag.total_mass;

			if ( std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) > 10000.0 ) {
				std::cout << "\n DIAGNOSTICS BENCHMARK DIFF H:\t" << "INF" << std::endl;
				//std::cout << "\n DIAGNOSTICS MASS DIFF:\t" << diagnostics_mass_start << " "<< simVars.diag.total_mass << " "<<std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
				std::cerr << "\n DIAGNOSTICS MASS DIFF TOO LARGE:\t" << std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
				exit(1);
			}

			// Print max abs difference of vars to initial conditions (this gives the error in steady state cases)
			//if ((simVars.setup.scenario >= 0 && simVars.setup.scenario <= 4) || simVars.setup.scenario == 13)
			//{
			// Height
			benchmark.t0_error_max_abs_h = (prog_h-t0_prog_h).reduce_maxAbs() ;
			o_ostream << "\t" << benchmark.t0_error_max_abs_h;

			// Velocity u
			benchmark.t0_error_max_abs_u = (prog_u-t0_prog_u).reduce_maxAbs();
			o_ostream << "\t" << benchmark.t0_error_max_abs_u;

			// Velocity v
			benchmark.t0_error_max_abs_v = (prog_v-t0_prog_v).reduce_maxAbs();
			o_ostream << "\t" << benchmark.t0_error_max_abs_v;
			//}

			if (param_compute_error && simVars.pde.use_nonlinear_equations==0)
			{
				compute_errors();

				o_ostream << "\t" << benchmark.analytical_error_rms_h;
				o_ostream << "\t" << benchmark.analytical_error_rms_u;
				o_ostream << "\t" << benchmark.analytical_error_rms_v;

				o_ostream << "\t" << benchmark.analytical_error_maxabs_h;
				o_ostream << "\t" << benchmark.analytical_error_maxabs_u;
				o_ostream << "\t" << benchmark.analytical_error_maxabs_v;
			}

			o_ostream << std::endl;
		}

		if (simVars.misc.output_each_sim_seconds > 0)
			while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;

		return true;
	}


public:
	void compute_errors()
	{
		// Compute exact solution for linear part and compare with numerical solution

		// Variables on unstaggered A-grid
		PlaneData t_h = t0_prog_h;
		PlaneData t_u = t0_prog_u;
		PlaneData t_v = t0_prog_v;

		if (simVars.disc.use_staggering) // Remap in case of C-grid
		{
			//remap initial condition to A grid
			//sampler2D.bilinear_scalar(t0_u, posx_a, posy_a, t_u, stag_u[0], stag_u[1]);
			//t_u=op.avg_f_x(t0_u); //equiv to bilinear
			sampler2D.bicubic_scalar(t0_prog_u, posa_x, posa_y, t_u, staggering.u[0], staggering.u[1]);

			//sampler2D.bilinear_scalar(t0_v, posx_a, posy_a, t_v, stag_v[0], stag_v[1]);
			//t_v=op.avg_f_y(t0_v); //equiv to bilinear
			sampler2D.bicubic_scalar(t0_prog_v, posa_x, posa_y, t_v, staggering.v[0], staggering.v[1]);
		}

		double o_dt;
		//Run exact solution for linear case
		timeSteppers.l_direct->run_timestep(
				t_h, t_u, t_v,
				o_dt,	// compute direct solution
				simVars.timecontrol.current_simulation_time,
				0,			// initial condition given at time 0
				simVars.timecontrol.max_simulation_time
		);

		// Analytical solution at specific time on orginal grid (stag or not)
		PlaneData ts_h = t_h;
		PlaneData ts_u = t_u;
		PlaneData ts_v = t_v;

		// Recover data in C grid using interpolations
		if (simVars.disc.use_staggering)
		{
			// Remap A grid to C grid

			//Temporary displacement for U points
			//sampler2D.bilinear_scalar(t_u, pos_x, tmp, ts_u, stag_h[0], stag_h[1]);
			//t_u=op.avg_b_x(t0_u); //equiv to bilinear
			sampler2D.bicubic_scalar(t_u, pos_x, pos_y, ts_u, staggering.h[0], staggering.h[1]+0.5);

			//Temporary displacement for V points
			//sampler2D.bilinear_scalar(t_v, tmp, pos_y, ts_v, stag_h[0], stag_h[1]);
			//t_v=op.avg_b_y(t0_v); //equiv to bilinear
			sampler2D.bicubic_scalar(t_v, pos_x, pos_y, ts_v, staggering.h[0]+0.5, staggering.h[1]);
		}

		benchmark.analytical_error_rms_h = (ts_h-prog_h).reduce_rms_quad();
		benchmark.analytical_error_rms_u = (ts_u-prog_u).reduce_rms_quad();
		benchmark.analytical_error_rms_v = (ts_v-prog_v).reduce_rms_quad();

		benchmark.analytical_error_maxabs_h = (ts_h-prog_h).reduce_maxAbs();
		benchmark.analytical_error_maxabs_u = (ts_u-prog_u).reduce_maxAbs();
		benchmark.analytical_error_maxabs_v = (ts_v-prog_v).reduce_maxAbs();

//		std::cout << std::endl;
//		std::cout << ts_h.reduce_rms_quad() << std::endl;

	}



public:
	bool should_quit()
	{
		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		if (simVars.timecontrol.max_simulation_time != -1 && simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time)
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
			{&prog_h,	"h"},
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
			PlaneData t_h = t0_prog_h;
			PlaneData t_u = t0_prog_u;
			PlaneData t_v = t0_prog_v;

			double o_dt;
			timeSteppers.l_direct->run_timestep(
					t_h, t_u, t_v,
					o_dt,
					simVars.timecontrol.current_simulation_time,
					0,
					simVars.timecontrol.max_simulation_time
				);

			switch(simVars.misc.vis_id)
			{
			case -1:
				vis = t_h;			//Exact solution
				break;

			case -2:
				vis = t_h-prog_h;	// difference to exact solution
				break;

			case -3:
				vis = t0_prog_h-prog_h;	// difference to initial condition
				break;
			}

			*o_dataArray = &vis;
			*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
			return;
		}

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
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
			prog_h.file_physical_saveData_ascii("swe_rexi_dump_h.csv");
			prog_u.file_physical_saveData_ascii("swe_rexi_dump_u.csv");
			prog_v.file_physical_saveData_ascii("swe_rexi_dump_v.csv");
			break;

		case 'C':
			// dump data arrays to VTK
			prog_h.file_physical_saveData_vtk("swe_rexi_dump_h.vtk", "Height");
			prog_u.file_physical_saveData_vtk("swe_rexi_dump_u.vtk", "U-Velocity");
			prog_v.file_physical_saveData_vtk("swe_rexi_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_h.file_physical_loadData("swe_rexi_dump_h.csv", simVars.setup.input_data_binary);
			prog_u.file_physical_loadData("swe_rexi_dump_u.csv", simVars.setup.input_data_binary);
			prog_v.file_physical_loadData("swe_rexi_dump_v.csv", simVars.setup.input_data_binary);
			break;
		}
	}
#endif


	bool instability_detected()
	{
		return !(	prog_h.reduce_boolean_all_finite() &&
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

		// use REXI
		timestepping_swe_plane_rexi.setup(
				0,
				simVars.rexi.rexi_h,
				simVars.rexi.rexi_M,
				simVars.rexi.rexi_L,

				simVars.disc.res_physical,
				simVars.sim.domain_size,
				simVars.rexi.rexi_use_half_poles,
				simVars.rexi.rexi_normalization
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


		*parareal_data_start.data_arrays[0] = prog_h;
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

		prog_h = *parareal_data_start.data_arrays[0];
		prog_u = *parareal_data_start.data_arrays[1];
		prog_v = *parareal_data_start.data_arrays[2];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		while (simVars.timecontrol.current_simulation_time != timeframe_end)
		{
			this->run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// copy to buffers
		*parareal_data_fine.data_arrays[0] = prog_h;
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

		prog_h = *parareal_data_start.data_arrays[0];
		prog_u = *parareal_data_start.data_arrays[1];
		prog_v = *parareal_data_start.data_arrays[2];

		timestepping_swe_plane_rexi.run_timestep_implicit_ts(
				prog_h, prog_u, prog_v,
				timeframe_end - timeframe_start,
				op,
				simVars
		);


		// copy to buffers
		*parareal_data_coarse.data_arrays[0] = prog_h;
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
			tmp = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			convergence = std::max(
					convergence,
					(*parareal_data_output.data_arrays[k]-tmp).reduce_maxAbs()
				);

			*parareal_data_output.data_arrays[k] = tmp;
		}

		simVars.timecontrol.current_simulation_time = timeframe_end;
		prog_h = *parareal_data_output.data_arrays[0];
		prog_u = *parareal_data_output.data_arrays[1];
		prog_v = *parareal_data_output.data_arrays[2];

		if (param_compute_error && simVars.pde.use_nonlinear_equations==0)
		{
			compute_errors();
			std::cout << "maxabs error compared to analytical solution: " << benchmark_analytical_error_maxabs_h << std::endl;
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
			FatalError()
		{
				std::cerr << "MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye..." << std::endl;
				exit(-1);
		}
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

	MemBlockAlloc::setup();

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"compute-error",
			"boundary-id",
			"initial-freq-x-mul",		/// frequency multipliers for special scenario setup
			"initial-freq-y-mul",
			"lin-exp-analyt",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 0;	// compute error - default no
	simVars.bogus.var[1] = 0; 	//boundary
	simVars.bogus.var[2] = 0;  //frequency in x for waves test case
	simVars.bogus.var[3] = 0;  //frequency in y for waves test case
	simVars.bogus.var[4] = 0;  // Use analytical linear operator exponential

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--timestepping-mode [0/1/2/3/4/5]	Timestepping method to use" << std::endl;
		std::cout << "	                            0: RKn with Finite-difference (default)" << std::endl;
		std::cout << "	                            1: REXI (SL-REXI if nonlinear)" << std::endl;
		std::cout << "	                            2: Direct solution in spectral space" << std::endl;
		std::cout << "	                            3: Implicit Euler Spectral " << std::endl;
		std::cout << "	                            4: Semi-Lagrangian Advection only" << std::endl;
		std::cout << "	                            5: Semi-Lagrangian Semi-implicit Spectral" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--compute-error [0/1]	Compute the errors" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--staggering [0/1]		Use staggered grid" << std::endl;
		std::cout << std::endl;
		std::cout << "	--boundary-id [0,1,...]	    Boundary id" << std::endl;
		std::cout << "                              0: no boundary (default)" << std::endl;
		std::cout << "                              1: centered box" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-zero-before-solving [0/1]	Zero the solution for the iterative solver (default=0)" << std::endl;
		std::cout << std::endl;
		//std::cout << "	--nonlinear [0/1/2]	   Form of equations:" << std::endl;
		//std::cout << "						     0: Linear SWE (default)" << std::endl;
		//std::cout << "						     1: Full nonlinear SWE" << std::endl;
		//std::cout << "						     2: Linear SWE + nonlinear advection only (needs -H to be set)" << std::endl;
		//std::cout << std::endl;
		std::cout << "	--lin-exp-analyt [0/1]	Use analytical exponential of linear operator (default=0)" << std::endl;
		std::cout << std::endl;


#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}
	// Calculate error flag
	param_compute_error = simVars.bogus.var[0];

	// Frequency for certain initial conditions
	param_initial_freq_x_mul = simVars.bogus.var[2];
	param_initial_freq_y_mul = simVars.bogus.var[3];

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	// Print header
	std::cout << std::endl;
	if (simVars.pde.use_nonlinear_equations == 1)
		std::cout << "Solving full nonlinear SW equations" << std::endl;
	else
		std::cout << "Solving linear SW equations" << std::endl;

	std::cout << "-----------------------------" << std::endl;
	simVars.outputConfig();
	std::cout << "Computing error: " << param_compute_error << std::endl;
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

			simulationSWE->reset();

			//Time counter
			Stopwatch time;

			//Diagnostic measures at initial stage
			//double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

			// Initialize diagnostics
			if (simVars.misc.verbosity > 0)
			{
				simulationSWE->update_diagnostics();
				diagnostics_energy_start = simVars.diag.total_energy;
				diagnostics_mass_start = simVars.diag.total_mass;
				diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
			}

#if SWEET_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			//Start counting time
			time.reset();


			if (simVars.disc.normal_mode_analysis_generation)
			{
				simulationSWE->normal_mode_analysis();
			}
			else
			{
				// Main time loop
				while(true)
				{
					if (simulationSWE->timestep_output(buf))
					{
						// string output data

						std::string output = buf.str();
						buf.str("");

						// This is an output printed on screen or buffered to files if > used
						std::cout << output;
					}

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
//			if (simVars.disc.timestepping_method != 0)
//				std::cout << "REXI alpha.size(): " << simulationSWE->timeStepperstimestepping_swe_plane_rexi.rexi.alpha.size() << std::endl;

			if (simVars.misc.verbosity > 0)
			{
				std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((simVars.diag.total_energy-diagnostics_energy_start)/diagnostics_energy_start) << std::endl;
				std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
				std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((simVars.diag.total_potential_enstrophy-diagnostics_potential_entrophy_start)/diagnostics_potential_entrophy_start) << std::endl;

				if (param_compute_error)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark.t0_error_max_abs_h << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark.t0_error_max_abs_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark.t0_error_max_abs_v << std::endl;
				}
			}

			if (param_compute_error && simVars.pde.use_nonlinear_equations == 0)
			{
				simulationSWE->compute_errors();
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
		if (simVars.disc.timestepping_method == SimulationVariables::Discretization::REXI)
		{
			PlaneOperators op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

			SWE_Plane_TS_l_rexi rexiSWE(simVars, op);

			/*
			 * Setup our little dog REXI
			 */
			rexiSWE.setup(
					simVars.rexi.rexi_h,
					simVars.rexi.rexi_M,
					simVars.rexi.rexi_L,
					simVars.rexi.rexi_use_half_poles,
					simVars.rexi.rexi_normalization
				);

			bool run = true;

			PlaneData prog_h(planeDataConfig);
			PlaneData prog_u(planeDataConfig);
			PlaneData prog_v(planeDataConfig);

			MPI_Barrier(MPI_COMM_WORLD);

			do
			{
				double o_dt;
				// REXI time stepping
				rexiSWE.run_timestep(
						prog_h,
						prog_u,
						prog_v,
						o_dt,					///< time step restriction
						-simVars.sim.CFL,		///< if this value is not equal to 0, use this time step size instead of computing one
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time

				);
			}
			while(!rexiSWE.final_timestep);
		}
	}
#endif


#if SWEET_MPI
	if (simVars.disc.timestepping_method > 0)
	{
		// synchronize REXI
		if (mpi_rank == 0)
			SWE_Plane_TS_l_rexi::MPI_quitWorkers(planeDataConfig);
	}

	MPI_Finalize();
#endif

	return 0;
}

