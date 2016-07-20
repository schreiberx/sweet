/*
 * SWE with nonlinear part using REXI test
 *
 */

#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/Sampler2D.hpp>
#include <sweet/SemiLagrangian.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "rexiswe/RexiSWE.hpp"



#ifndef SWEET_MPI
#	define SWEET_MPI 1
#endif

#if SWEET_MPI
#	include <mpi.h>
#endif

#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif
#include <parareal/Parareal.hpp>


// Input parameters (cmd line)

//general parameters
SimulationVariables simVars;

//specific parameters
double param_rexi_h;
double param_rexi_m;
double param_rexi_l;
bool param_rexi_half;
int param_timestepping_mode;
bool param_compute_error;
bool param_use_staggering;
bool param_rexi_use_spectral_differences_for_complex_array;
int param_rexi_helmholtz_solver_id;
double param_rexi_helmholtz_solver_eps;
bool param_rexi_zero_before_solving;
int param_boundary_id;
int param_nonlinear;

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
	DataArray<2> prog_h, prog_u, prog_v;

	// beta plane
	DataArray<2> beta_plane;

	// Prognostic variables at time step t-dt
	DataArray<2> prog_h_prev, prog_u_prev, prog_v_prev;

	// Diagnostics - Vorticity and potential vorticity
	DataArray<2> eta, q;

	//visualization variable
	DataArray<2> vis;

	// temporary variables - may be overwritten, use locally
	DataArray<2> tmp, tmp0, tmp1, tmp2, tmp3;

	// Variables to keep track of boundary
	DataArray<2> boundary_mask;
	DataArray<2> boundary_mask_inv;

	// Initial values for comparison with analytical solution
	DataArray<2> t0_prog_h, t0_prog_u, t0_prog_v;

	// Forcings
	DataArray<2> force_h, force_u, force_v;

	// Nonlinear terms relative to h, u, v and its previous values
	//    These have already exp(dtL/2) applied to them.
	DataArray<2> N_h, N_u, N_v;
	DataArray<2> N_h_prev, N_u_prev, N_v_prev;

	// Points mapping [0,simVars.sim.domain_size[0])x[0,simVars.sim.domain_size[1])
	// with resolution simVars.sim.resolution
	DataArray<2> pos_x, pos_y;

	// Arrival points for semi-lag
	DataArray<2> posx_a, posy_a;

	// Departure points for semi-lag
	DataArray<2> posx_d, posy_d;


	//Staggering displacement array (use 0.5 for each displacement)
	// [0] - delta x of u variable
	// [1] - delta y of u variable
	// [2] - delta x of v variable
	// [3] - delta y of v variable
	// Default - A grid (there is shift in x,y of 1/2 to all vars)
	// For C grid use {0,-0.5,-0.5,0}
	double stag_displacement[4] = {-0.5,-0.5,-0.5,-0.5};
	double stag_h[2] = {-0.5,-0.5};
	double stag_u[2] = {-0.5,-0.5};
	double stag_v[2] = {-0.5,-0.5};

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	//Max difference to initial conditions
	double benchmark_diff_h;
	double benchmark_diff_u;
	double benchmark_diff_v;

	// Error measures L2 norm
	double benchmark_analytical_error_rms_h;
	double benchmark_analytical_error_rms_u;
	double benchmark_analytical_error_rms_v;

	// Error measures max norm
	double benchmark_analytical_error_maxabs_h;
	double benchmark_analytical_error_maxabs_u;
	double benchmark_analytical_error_maxabs_v;

	// Finite difference operators
	Operators2D op;

	// Runge-Kutta stuff
	TimesteppingRK timestepping;

	// Rexi stuff
	RexiSWE rexiSWE;

	// Interpolation stuff
	Sampler2D sampler2D;

	// Semi-Lag stuff
	SemiLagrangian semiLagrangian;

public:
	SimulationInstance()	:
	// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_h(simVars.disc.res),
		prog_u(simVars.disc.res),
		prog_v(simVars.disc.res),

		beta_plane(simVars.disc.res),

		prog_h_prev(simVars.disc.res),
		prog_u_prev(simVars.disc.res),
		prog_v_prev(simVars.disc.res),

		eta(simVars.disc.res),
		q(simVars.disc.res),

		vis(simVars.disc.res),

		tmp(simVars.disc.res),
		tmp0(simVars.disc.res),
		tmp1(simVars.disc.res),
		tmp2(simVars.disc.res),
		tmp3(simVars.disc.res),

		boundary_mask(simVars.disc.res),
		boundary_mask_inv(simVars.disc.res),

		t0_prog_h(simVars.disc.res),
		t0_prog_u(simVars.disc.res),
		t0_prog_v(simVars.disc.res),

		force_h(simVars.disc.res),
		force_u(simVars.disc.res),
		force_v(simVars.disc.res),

		// @Martin: This should be added only when nonlinear model is ran. How to do that, since I can't put an if in the constructor?
		N_h(simVars.disc.res),
		N_u(simVars.disc.res),
		N_v(simVars.disc.res),

		N_h_prev(simVars.disc.res),
		N_u_prev(simVars.disc.res),
		N_v_prev(simVars.disc.res),

		pos_x(simVars.disc.res),
		pos_y(simVars.disc.res),

		posx_a(simVars.disc.res),
		posy_a(simVars.disc.res),

		posx_d(simVars.disc.res),
		posy_d(simVars.disc.res),

		//h_t(simVars.disc.res),

		// Initialises operators
		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
#if SWEET_PARAREAL != 0
		,
		_parareal_data_start_h(simVars.disc.res), _parareal_data_start_u(simVars.disc.res), _parareal_data_start_v(simVars.disc.res),
		_parareal_data_fine_h(simVars.disc.res), _parareal_data_fine_u(simVars.disc.res), _parareal_data_fine_v(simVars.disc.res),
		_parareal_data_coarse_h(simVars.disc.res), _parareal_data_coarse_u(simVars.disc.res), _parareal_data_coarse_v(simVars.disc.res),
		_parareal_data_output_h(simVars.disc.res), _parareal_data_output_u(simVars.disc.res), _parareal_data_output_v(simVars.disc.res),
		_parareal_data_error_h(simVars.disc.res), _parareal_data_error_u(simVars.disc.res), _parareal_data_error_v(simVars.disc.res)
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
		// Initialise diagnostics
		last_timestep_nr_update_diagnostics = -1;

		benchmark_diff_h = 0;
		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		benchmark_analytical_error_rms_h = 0;
		benchmark_analytical_error_rms_u = 0;
		benchmark_analytical_error_rms_v = 0;

		benchmark_analytical_error_maxabs_h = 0;
		benchmark_analytical_error_maxabs_u = 0;
		benchmark_analytical_error_maxabs_v = 0;

		simVars.timecontrol.current_timestep_nr = 0;
		simVars.timecontrol.current_simulation_time = 0;

		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_SPECTRAL_SPACE
		prog_h.set_spec_all(0, 0);
		prog_u.set_spec_all(0, 0);
		prog_v.set_spec_all(0, 0);
		boundary_mask.set_spec_all(0, 0);
#endif

		//Setup prog vars
		prog_h.set_all(simVars.setup.h0);
		prog_u.set_all(0);
		prog_v.set_all(0);
		boundary_mask.set_all(0);

		//Check if input parameters are adequate for this simulation
		if (param_use_staggering && simVars.disc.use_spectral_basis_diffs)
		{
			std::cerr << "Staggering and spectral basis not supported!" << std::endl;
			exit(1);
		}

		if (param_use_staggering && ( param_timestepping_mode != 0 && param_timestepping_mode != 4 ) )
		{
			std::cerr << "Staggering only supported for standard time stepping mode 0 and 4!" << std::endl;
			exit(1);
		}
		if (param_use_staggering && param_compute_error)
			std::cerr << "Warning: Staggered data will be interpolated to/from A-grid for exact linear solution" << std::endl;

		if(param_nonlinear > 0 && param_compute_error)
			std::cout << "Warning: Exact solution not possible in general for nonlinear swe. Using exact solution for linear case instead." << std::endl;

		if (param_use_staggering)
		{
			/*
			 *              ^
			 *              |
			 *       ______v0,1_____
			 *       |             |
			 *       |			   |
			 *       |   (0.5,0.5) |
			 *  u0,0 |->  H/P0,0   |u1,0 ->
			 *(0,0.5)|	           |
			 *       |      ^      |
			 *   q0,0|______|______|
			 * (0,0)      v0,0
			 *           (0.5,0)
			 *
			 *
			 * These staggering should be used when interpolating from a staggered variable
			 * If interpolating from A grid to C staggered, use negative of displacements.
			 *
			 */
			stag_displacement[0] = 0.0;  // u_dx
			stag_displacement[1] = -0.5; // u_dy
			stag_displacement[2] = -0.5; // v_dx
			stag_displacement[3] = 0.0;  // v_dy
			stag_h[0] = -0.5;
			stag_h[1] = -0.5;
			stag_u[0] = -0.0;
			stag_u[1] = -0.5;
			stag_v[0] = -0.5;
			stag_v[1] = -0.0;
		}

		//Setup sampler for future interpolations
		sampler2D.setup(simVars.sim.domain_size, simVars.disc.res);

		//Setup semi-lag
		semiLagrangian.setup(simVars.sim.domain_size, simVars.disc.res);

		//Setup general (x,y) grid with position points
		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
		    	/* Equivalent to q position on C-grid */
				pos_x.set(j, i, ((double)i)*simVars.sim.domain_size[0]/simVars.disc.res[0]); //*simVars.sim.domain_size[0];
				pos_y.set(j, i, ((double)j)*simVars.sim.domain_size[1]/simVars.disc.res[1]); //*simVars.sim.domain_size[1];
				//std::cout << i << " " << j << " " << pos_x.get(j,i) << std::endl;
			}
		}

		//Initialize arrival points with h position
		posx_a = pos_x+0.5*simVars.disc.cell_size[0];
		posy_a = pos_y+0.5*simVars.disc.cell_size[1];

		//std::cout << std::endl;
		//std::cout << "posx_a: " << posx_a.array_data_cartesian_space_valid << std::endl;
		//std::cout << std::endl;

		//@martin : Shouldn't this be in SWEValidationBenchmarks?
		auto return_h = [] (
				SimulationVariables &i_parameters,
				double x,
				double y
		) -> double
		{
			if (param_initial_freq_x_mul == 0)
				return SWEValidationBenchmarks::return_h(simVars, x, y);

			// Waves scenario
			// Remember to set up initial_freq_x_mul and initial_freq_y_mul
			double dx = x/i_parameters.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
			double dy = y/i_parameters.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
			return std::sin(2.0*dx)*std::cos(2.0*dy) - (1.0/5.0)*std::cos(2.0*dx)*std::sin(4.0*dy) + i_parameters.setup.h0;
		};


		auto return_u = [] (
				SimulationVariables &i_parameters,
				double x,
				double y
		) -> double
		{
			if (param_initial_freq_x_mul == 0)
				return SWEValidationBenchmarks::return_u(simVars, x, y);

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
				return SWEValidationBenchmarks::return_v(simVars, x, y);

			double dx = x/i_parameters.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
			double dy = y/i_parameters.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
			return std::cos(2.0*dx)*std::cos(4.0*dy);
		};

		// Set initial conditions given from SWEValidationBenchmarks
		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
				if (param_use_staggering) // C-grid
				{
					{
						// h
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						prog_h.set(j, i, return_h(simVars, x, y));
						t0_prog_h.set(j, i, return_h(simVars, x, y));
						force_h.set(j, i, SWEValidationBenchmarks::return_force_h(simVars, x, y));

						std::cerr << "WARNING: BETA PLANE ON C-GRID NOT SUPPORTED!" << std::endl;
						beta_plane.set(j, i, SWEValidationBenchmarks::return_f(simVars, x, y));
					}

					{
						// u space
						double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						prog_u.set(j,i, return_u(simVars, x, y));
						t0_prog_u.set(j, i, return_u(simVars, x, y));
						force_u.set(j, i, SWEValidationBenchmarks::return_force_u(simVars, x, y));
					}

					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						prog_v.set(j, i, return_v(simVars, x, y));
						t0_prog_v.set(j, i, return_v(simVars, x, y));
						force_v.set(j, i, SWEValidationBenchmarks::return_force_v(simVars, x, y));
					}
				}
				else // A-Grid (colocated grid)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

					prog_h.set(j, i, return_h(simVars, x, y));
					prog_u.set(j, i, return_u(simVars, x, y));
					prog_v.set(j, i, return_v(simVars, x, y));

					// beta plane
					if (simVars.sim.beta < 0)
					{
						double y_beta = (((double)j+0.5)/(double)simVars.disc.res[1]);
						beta_plane.set(j, i, 0);
					}
					else
					{
						double y_beta = (((double)j+0.5)/(double)simVars.disc.res[1]);
						beta_plane.set(j, i, simVars.sim.f0+simVars.sim.beta*y_beta);
					}

					t0_prog_h.set(j, i, return_h(simVars, x, y));
					t0_prog_u.set(j, i, return_u(simVars, x, y));
					t0_prog_v.set(j, i, return_v(simVars, x, y));

					force_h.set(j, i, SWEValidationBenchmarks::return_force_h(simVars, x, y));
					force_u.set(j, i, SWEValidationBenchmarks::return_force_u(simVars, x, y));
					force_v.set(j, i, SWEValidationBenchmarks::return_force_v(simVars, x, y));
				}
			}
		}

		//Initialise t-dt time step with initial condition
		prog_h_prev = prog_h;
		prog_u_prev = prog_u;
		prog_v_prev = prog_v;

		//Nonlinear variables
		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_SPECTRAL_SPACE
		N_h.set_spec_all(0, 0);
		N_u.set_spec_all(0, 0);
		N_v.set_spec_all(0, 0);
#endif
		N_h.set_all(0);
		N_u.set_all(0);
		N_v.set_all(0);

		N_h_prev=N_h;
		N_u_prev=N_u;
		N_v_prev=N_v;

		// Set boundary stuff
		if (param_boundary_id != 0)
		{
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					int boundary_flag = 0;
					switch(param_boundary_id)
					{
					case 0:
						break;

					case 1:
						boundary_flag =
								(i >= simVars.disc.res[0]*1/4) &&
								(i < simVars.disc.res[0]*3/4) &&
								(j >= simVars.disc.res[1]*1/4) &&
								(j < simVars.disc.res[1]*3/4)
							;
						break;


					case 2:
						boundary_flag =
								(i >= simVars.disc.res[0]*3/8) &&
								(i < simVars.disc.res[0]*5/8) &&
								(j >= simVars.disc.res[1]*3/8) &&
								(j < simVars.disc.res[1]*5/8)
							;
						break;


					case 3:
						boundary_flag |=
								(i >= simVars.disc.res[0]*1/8) &&
								(i < simVars.disc.res[0]*3/8) &&
								(j >= simVars.disc.res[1]*1/8) &&
								(j < simVars.disc.res[1]*3/8)
							;
						boundary_flag |=
								(i >= simVars.disc.res[0]*1/8) &&
								(i < simVars.disc.res[0]*3/8) &&
								(j >= simVars.disc.res[1]*5/8) &&
								(j < simVars.disc.res[1]*7/8)
							;
						boundary_flag |=
								(i >= simVars.disc.res[0]*5/8) &&
								(i < simVars.disc.res[0]*7/8) &&
								(j >= simVars.disc.res[1]*5/8) &&
								(j < simVars.disc.res[1]*7/8)
							;
						boundary_flag |=
								(i >= simVars.disc.res[0]*5/8) &&
								(i < simVars.disc.res[0]*7/8) &&
								(j >= simVars.disc.res[1]*	1/8) &&
								(j < simVars.disc.res[1]*3/8)
							;
						break;



					case 4:
						boundary_flag |=
								(i >= simVars.disc.res[0]*0/8) &&
								(i < simVars.disc.res[0]*2/8) &&
								(j >= simVars.disc.res[1]*0/8) &&
								(j < simVars.disc.res[1]*2/8)
							;
						boundary_flag |=
								(i >= simVars.disc.res[0]*1/8) &&
								(i < simVars.disc.res[0]*3/8) &&
								(j >= simVars.disc.res[1]*5/8) &&
								(j < simVars.disc.res[1]*7/8)
							;
						boundary_flag |=
								(i >= simVars.disc.res[0]*5/8) &&
								(i < simVars.disc.res[0]*7/8) &&
								(j >= simVars.disc.res[1]*6/8) &&
								(j < simVars.disc.res[1]*8/8)
							;
						boundary_flag |=
								(i >= simVars.disc.res[0]*5/8) &&
								(i < simVars.disc.res[0]*7/8) &&
								(j >= simVars.disc.res[1]*	1/8) &&
								(j < simVars.disc.res[1]*3/8)
							;
						break;

					default:
						std::cerr << "Unknown boundary id " << param_boundary_id << std::endl;
						exit(1);
					}

					boundary_mask.set(j, i, boundary_flag);
					boundary_mask_inv.set(j, i, 1.0-boundary_flag);
				}
			}
		}

		// Load data, if requested
		if (simVars.setup.input_data_filenames.size() > 0)
			prog_h.file_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 1)
		{
			prog_u.file_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);
			if (param_use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

		if (simVars.setup.input_data_filenames.size() > 2)
		{
			prog_v.file_loadData(simVars.setup.input_data_filenames[2].c_str(), simVars.setup.input_data_binary);
			if (param_use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}


		// Print info for REXI and setup REXI
		if (param_timestepping_mode == 1) // This is not necessary: || param_timestepping_mode == 3)
		{
			if (simVars.misc.verbosity > 0)
			{
				std::cout << "REXI: Using REXI for time integration" << std::endl;
				std::cout << "REXI: Using M=" << param_rexi_m << ", L=" << param_rexi_l << ", poles and sampling size h=" << param_rexi_h << std::endl;
			}

			if (simVars.sim.CFL >= 0)
			{
				std::cout << "Only constant time step size supported with REXI, use negative CFL to set constant time step size" << std::endl;
				exit(1);
			}

			// use REXI
			rexiSWE.setup(
					param_rexi_h,
					param_rexi_m,
					param_rexi_l,
					simVars.disc.res,
					simVars.sim.domain_size,
					param_rexi_half,
					param_rexi_use_spectral_differences_for_complex_array,
					param_rexi_helmholtz_solver_id,
					param_rexi_helmholtz_solver_eps
				);

			if (simVars.misc.verbosity > 2)
			{
				std::cout << "ALPHA:" << std::endl;
				for (std::size_t n = 0; n < rexiSWE.rexi.alpha.size(); n++)
					std::cout << rexiSWE.rexi.alpha[n] << std::endl;

				std::cout << "BETA:" << std::endl;
				for (std::size_t n = 0; n < rexiSWE.rexi.beta_re.size(); n++)
					std::cout << rexiSWE.rexi.beta_re[n] << std::endl;
			}
		}

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

		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res[0]*(double)simVars.disc.res[1]);


		// mass
		simVars.diag.total_mass = prog_h.reduce_sum_quad() * normalization;

		// energy
		simVars.diag.total_energy = 0.5*((
				prog_h*prog_h +
				prog_h*prog_u*prog_u +
				prog_h*prog_v*prog_v
			).reduce_sum_quad()) * normalization;

		// potential verticity and pot. enstropy
		eta = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + simVars.sim.f0) / prog_h;
		simVars.diag.total_potential_enstrophy = 0.5*(eta*eta*prog_h).reduce_sum_quad() * normalization;
	}


	void compute_upwinding_P_updates(
			const DataArray<2> &i_h,		///< prognostic variables (at T=tn)
			const DataArray<2> &i_u,		///< prognostic variables (at T=tn+dt)
			const DataArray<2> &i_v,		///< prognostic variables (at T=tn+dt)

			DataArray<2> &o_P_t				///< time updates (at T=tn+dt)
	)
	{
		std::cerr << "TODO: implement, is this really possible for non-staggered grid? (averaging of velocities required)" << std::endl;
		exit(-1);
		//             |                       |                       |
		// --v---------|-----------v-----------|-----------v-----------|
		//   h-1       u0          h0          u1          h1          u2
		//

		// same as above, but formulated in a finite-difference style
		o_P_t =
			(
				(
					// u is positive
					simVars.setup.h0*i_u.return_value_if_positive()	// inflow
					-simVars.setup.h0*op.shift_left(i_u.return_value_if_positive())					// outflow

					// u is negative
					+(simVars.setup.h0*i_u.return_value_if_negative())	// outflow
					-op.shift_left(simVars.setup.h0*i_u.return_value_if_negative())		// inflow
				)*(1.0/simVars.disc.cell_size[0])	// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					simVars.setup.h0*i_v.return_value_if_positive()		// inflow
					-simVars.setup.h0*op.shift_down(i_v.return_value_if_positive())					// outflow

					// v is negative
					+(simVars.setup.h0*i_v.return_value_if_negative())	// outflow
					-op.shift_down(simVars.setup.h0*i_v.return_value_if_negative())	// inflow
				)*(1.0/simVars.disc.cell_size[1])
			);
	}



	void boundary_action()
	{

		if (param_boundary_id == 0)
			return;

		if (param_use_staggering)
		{
			/*
			 * 0     1       1       0       0
			 *       h0      h1      h2      h3
			 *   |---+---|---+---|---+---|---+---|
			 *       |XXXXXXX|
			 *
			 *   u0      u1      u2      u3      u4
			 *
			 *   0       1       1       0       0	// boundary mask at u
			 *   1       0       0       1       1	// inv boundary mask at u
			 */
			// process U
			prog_u =
					// velocities outside of boundary
					prog_u*(boundary_mask_inv*op.shift_right(boundary_mask_inv))
					// left boundary
					- op.shift_right(prog_u*boundary_mask_inv)*boundary_mask
					// right boundary
					- op.shift_left(prog_u)*boundary_mask_inv*op.shift_right(boundary_mask);

			prog_v =
					// velocities outside of boundary
					prog_v*(boundary_mask_inv*op.shift_up(boundary_mask_inv))
					// left boundary
					- op.shift_up(prog_v*boundary_mask_inv)*boundary_mask
					// right boundary
					- op.shift_down(prog_v)*boundary_mask_inv*op.shift_up(boundary_mask);

			prog_h = boundary_mask_inv*prog_h + boundary_mask*simVars.setup.h0;
		}
		else
		{
			//          |XXXXXXXXXXXXXXX
			//  u0     u1    u2     u3
			//  0      1     1      1	// boundary
			//  1      0     0      0	// inv
#if 1
			prog_u = boundary_mask_inv*prog_u
					+ op.shift_right(op.shift_right(boundary_mask_inv*prog_u)*boundary_mask)
					+ op.shift_right(boundary_mask_inv*prog_u)*boundary_mask
					+ op.shift_left(op.shift_left(boundary_mask_inv*prog_u)*boundary_mask)
					+ op.shift_left(boundary_mask_inv*prog_u)*boundary_mask
					;
#endif
#if 1
			prog_v = boundary_mask_inv*prog_v
					+ op.shift_up(op.shift_up(boundary_mask_inv*prog_v)*boundary_mask)
					+ op.shift_up(boundary_mask_inv*prog_v)*boundary_mask
					+ op.shift_down(op.shift_down(boundary_mask_inv*prog_v)*boundary_mask)
					+ op.shift_down(boundary_mask_inv*prog_v)*boundary_mask
					;
#endif

			prog_h = boundary_mask_inv*prog_h + boundary_mask*simVars.setup.h0;
		}
	}

	// Main routine for method to be used in case of finite differences
	void p_run_euler_timestep_update(
			const DataArray<2> &i_h,	///< prognostic variables
			const DataArray<2> &i_u,	///< prognostic variables
			const DataArray<2> &i_v,	///< prognostic variables

			DataArray<2> &o_h_t,	///< time updates
			DataArray<2> &o_u_t,	///< time updates
			DataArray<2> &o_v_t,	///< time updates

			double &o_dt,			///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{

		/*
		 * TIME STEP SIZE
		 */
		if (i_fixed_dt > 0)
		{
			o_dt = i_fixed_dt;
		}
		else
		{
			/*
			 * If the timestep size parameter is negative, we use the absolute value of this one as the time step size
			 */
			if (i_fixed_dt < 0)
			{
				o_dt = -i_fixed_dt;
			}
			else
			{
#if 0
				double limit_speed = std::min(simVars.disc.cell_size[0]/i_u.reduce_maxAbs(), simVars.disc.cell_size[1]/i_v.reduce_maxAbs())*(1.0/simVars.sim.f0);

				if (std::abs(simVars.sim.f0) == 0)
					limit_speed = std::numeric_limits<double>::infinity();

				//	double hx = simVars.disc.cell_size[0];
				//	double hy = simVars.disc.cell_size[1];

				// limit by viscosity
				double limit_visc = std::numeric_limits<double>::infinity();
				/*
				if (simVars.sim.viscosity > 0)
					limit_visc = (hx*hx*hy*hy)/(4.0*simVars.sim.viscosity*simVars.sim.viscosity);
				if (simVars.sim.viscosity_order > 0)
					limit_visc = std::min((hx*hx*hx*hx*hy*hy*hy*hy)/(16.0*simVars.sim.viscosity_order*simVars.sim.viscosity_order), limit_visc);
				 */
				// limit by gravitational acceleration
				double limit_gh = std::min(simVars.disc.cell_size[0], simVars.disc.cell_size[1])/std::sqrt(simVars.sim.g*simVars.setup.h0);

				if (simVars.misc.verbosity > 2)
					std::cerr << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

				o_dt = simVars.sim.CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);
#endif

				double k = 1.0/std::min(simVars.disc.cell_size[0], simVars.disc.cell_size[1]);
				o_dt = simVars.sim.CFL /
					std::sqrt(
						(k*k)*simVars.sim.g*simVars.setup.h0
						+
						simVars.sim.f0*simVars.sim.f0
					);
			}
		}

		// A- grid method
		if (!param_use_staggering)
		{
			if (param_nonlinear > 0)	 // nonlinear case
			{
				std::cout << "Only linear swe are setup for unstaggered grids " << std::endl;
				exit(1);
			}

			boundary_action();

			/*
			 * linearized non-conservative (advective) formulation:
			 *
			 * h_t = -h0*u_x - h0*v_ym
			 * u_t = -g * h_x + f*v
			 * v_t = -g * h_y - f*u
			 */

			o_u_t = -simVars.sim.g*op.diff_c_x(i_h);
			o_v_t = -simVars.sim.g*op.diff_c_y(i_h);
			if (simVars.sim.beta == 0.0)
			{
				o_u_t += simVars.sim.f0*i_v;
				o_v_t -= simVars.sim.f0*i_u;
			}
			else
			{
				o_u_t += beta_plane*i_v;
				o_v_t -= beta_plane*i_u;
			}


			if (simVars.sim.viscosity != 0)
			{
				o_u_t -= op.diffN_x(i_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
				o_v_t -= op.diffN_y(i_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			}

			boundary_action();

			if (!simVars.disc.timestepping_up_and_downwinding)
			{
				// standard update
				o_h_t = -(op.diff_c_x(i_u) + op.diff_c_y(i_v))*simVars.setup.h0;
			}
			else
			{
				// up/down winding
				compute_upwinding_P_updates(
						i_h,
						i_u,
						i_v,
						o_h_t
					);
			}

			boundary_action();

#if 0
			if (simVars.sim.potential_viscosity != 0)
				o_h_t -= op.diff2(i_h)*simVars.sim.potential_viscosity;

			if (simVars.sim.potential_hyper_viscosity != 0)
				o_h_t -= op.diff4(i_h)*simVars.sim.potential_hyper_viscosity;
#endif
		}
		else // param_use_staggering = true
		{
			boundary_action();

			// STAGGERED GRID

			DataArray<2>& U = tmp0;
			DataArray<2>& V = tmp1;
			DataArray<2>& H = tmp2;

			/* Sadourny energy conserving scheme
			 *
			 * Note, that this grid does not follow the formulation
			 * in the paper of Robert Sadourny, but looks as follows:
			 *
			 *              ^
			 *              |
			 *       ______v0,1_____
			 *       |             |
			 *       |			   |
			 *       |             |
			 *  u0,0 |->  H/P0,0   |u1,0 ->
			 *(0,0.5)|			   |
			 *       |      ^      |
			 *   q0,0|______|______|
			 * (0,0)      v0,0
			 *           (0.5,0)
			 *
			 * V_t + q N x (P V) + grad( g P + 1/2 V*V) = 0
			 * P_t + div(P V) = 0
			 */
			/*
			 * U and V updates
			 */
			if(param_nonlinear > 0) //nonlinear case
			{
				U = op.avg_b_x(i_h)*i_u;
				V = op.avg_b_y(i_h)*i_v;
			}
			else // linear case
			{
				U = simVars.setup.h0*i_u;
				V = simVars.setup.h0*i_v;
			}

			if(param_nonlinear > 0) //nonlinear case
				H = simVars.sim.g*i_h + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));
			else //linear case
				H = simVars.sim.g*i_h;// + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));


			if(param_nonlinear > 0) //nonlinear case
			{
				// Potential vorticity
				if(op.avg_b_x(op.avg_b_y(i_h)).reduce_min() < 0.00000001)
				{
					std::cerr << "Test case not adequate for vector invariant formulation. Null or negative water height" << std::endl;
					exit(1);
				}
				q = (op.diff_b_x(i_v) - op.diff_b_y(i_u) + simVars.sim.f0) / op.avg_b_x(op.avg_b_y(i_h));

				// u, v tendencies
				// Energy conserving scheme
				o_u_t = op.avg_f_y(q*op.avg_b_x(V)) - op.diff_b_x(H);
				o_v_t = -op.avg_f_x(q*op.avg_b_y(U)) - op.diff_b_y(H);
			}
			else //linear case
			{
				o_u_t = op.avg_f_y(simVars.sim.f0*op.avg_b_x(i_v)) - op.diff_b_x(H);
				o_v_t = -op.avg_f_x(simVars.sim.f0*op.avg_b_y(i_u)) - op.diff_b_y(H);
			}


			/*
			 * VISCOSITY
			 */
			if (simVars.sim.viscosity != 0)
			{
				o_u_t -= op.diffN_x(i_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
				o_v_t -= op.diffN_y(i_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			}

			/*
			 * P UPDATE
			 */
			if (!simVars.disc.timestepping_up_and_downwinding)
			{
				if(param_nonlinear > 0){ //full nonlinear divergence
					// standard update
					o_h_t = -op.diff_f_x(U) - op.diff_f_y(V);
				}
				else //use linear divergence
				{
					o_h_t = -op.diff_f_x(simVars.setup.h0*i_u) - op.diff_f_y(simVars.setup.h0*i_v);
				}
			}
			else
			{
				// up/down winding
				compute_upwinding_P_updates(
						i_h,
						i_u,
						i_v,
						o_h_t
					);
			}
		}
	}



#if 0
	// doesn't seem to be necessary
	/**
	 * wrapper to unify interfaces for REXI and analytical solution for L operator
	 */
	void rexi_run_timestep_wrapper(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		double i_local_timestep_size
	)
	{
		switch(param_timestepping_mode)
		{
		case 1:
			// REXI time stepping
			assert(simVars.sim.CFL < 0);
			rexiSWE.run_timestep(
					prog_h, prog_u, prog_v,
					i_local_timestep_size,
					op,
					simVars,
					param_rexi_zero_before_solving
			);
			break;

		case 2:
			if (param_use_staggering)
			{
				std::cerr << "Direct solution on staggered grid not supported!" << std::endl;
				exit(1);
			}

			// Analytical solution
			rexiSWE.run_timestep_direct_solution(
					prog_h, prog_u, prog_v,
					i_local_timestep_size,
					op,
					simVars
			);
			break;

		default:
			assert(false);
			std::cerr << "Timestepping method in this wrapper not supported" << std::endl;
			exit(-1);
		}
	}
#endif



	/**
	 * Execute a single simulation time step
	 */
	void run_timestep(
//			double i_max_simulation_time
	)
	{
		double o_dt;

		if (param_timestepping_mode == 0)
		{
			// either set time step size to 0 for autodetection or to
			// a positive value to use a fixed time step size
			simVars.timecontrol.current_timestep_size = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 0);

			// standard time stepping
			timestepping.run_rk_timestep(
					this,
					&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					prog_h, prog_u, prog_v,
					o_dt,
					simVars.timecontrol.current_timestep_size,
					simVars.disc.timestepping_runge_kutta_order,
					simVars.timecontrol.current_simulation_time,
					simVars.timecontrol.max_simulation_time
				);
		}
		else if (param_timestepping_mode == 1) //REXI
		{
			assert(simVars.sim.CFL < 0);
			o_dt = -simVars.sim.CFL;

			// REXI time stepping for nonlinear eq - semi-lagrangian scheme (SL-REXI)
			if (param_nonlinear>0)
			{

				//Calculate departure points
				semiLagrangian.semi_lag_departure_points_settls(
								prog_u_prev, prog_v_prev,
								prog_u,	prog_v,
								posx_a,	posy_a,
								o_dt,
								posx_d,	posy_d,
								stag_displacement
						);

				//std::cout<<"Departure points movement (x,y): " << (posx_a-posx_d).reduce_maxAbs() << " " <<(posy_a-posy_d).reduce_maxAbs() << std::endl;

				//Save current fields for next time step
				prog_u_prev = prog_u;
				prog_v_prev = prog_v;
				prog_h_prev = prog_h;


				DataArray<2>& U = tmp0;
				DataArray<2>& V = tmp1;
				DataArray<2>& H = tmp2;

				U = prog_u;
				V = prog_v;
				H = prog_h;

				/*
				 *
				std::cout<<"U"<<std::endl;
				U.printArrayData();
				std::cout<<"V"<<std::endl;
				V.printArrayData();
				std::cout<<"H"<<std::endl;
				H.printArrayData();

				std::cout<<"prog_u"<<std::endl;
				prog_u.printArrayData();
				std::cout<<"prog_v"<<std::endl;
				prog_v.printArrayData();
				std::cout<<"prog_h"<<std::endl;
				prog_h.printArrayData();
				std::cout<<std::endl;
				*/

				if(param_nonlinear==1) //Full nonlinear swe
				{
					// First calculate the linear part
					rexiSWE.run_timestep( H, U, V, o_dt, op, simVars, param_rexi_zero_before_solving);
					//rexiSWE.run_timestep_direct_solution( H, U, V, o_dt, op, simVars );

					//Zero nonlinear vectors
#if SWEET_USE_SPECTRAL_SPACE && 0
							N_h.set_spec_all(0, 0);
							N_u.set_spec_all(0, 0);
							N_v.set_spec_all(0, 0);
#endif
							N_h.set_all(0);
							N_u.set_all(0);
							N_v.set_all(0);

					//Calculate divergence spectrally or with finite differences
					//if(simVars.disc.use_spectral_basis_diffs) //Spectral derivatives
					//	N_h=-prog_h*(op.diff_c_x(prog_u) + op.diff_c_y(prog_v));
					//else //Finite differences - needs further averaging for A grid?
					N_h=-prog_h*(op.diff_c_x(prog_u) + op.diff_c_y(prog_v));

					//Calculate exp(Ldt/2 N(u))
					rexiSWE.run_timestep( N_h, N_u, N_v, o_dt/2.0, op, simVars, param_rexi_zero_before_solving);

					//Use previous step to calculate main term to be interpolated
					H=H+(o_dt)*N_h-(o_dt/2.0)*N_h_prev;
					U=U+(o_dt)*N_u-(o_dt/2.0)*N_u_prev;
					V=V+(o_dt)*N_v-(o_dt/2.0)*N_v_prev;

					//Now interpolate to the the departure points
					//Departure points are set for physical space
					sampler2D.bicubic_scalar( H, posx_d, posy_d, prog_h, stag_h[0],	stag_h[1]);
					sampler2D.bicubic_scalar( U, posx_d, posy_d, prog_u, stag_u[0], stag_u[1]);
					sampler2D.bicubic_scalar( V, posx_d, posy_d, prog_v, stag_v[0], stag_v[1]);

					//Add nonlinear part attributed to arrival points
					prog_h=prog_h+(o_dt/2.0)*N_h;
					prog_u=prog_u+(o_dt/2.0)*N_u;
					prog_v=prog_v+(o_dt/2.0)*N_v;

					//Save current nonlinear values for next time step
					N_h_prev=N_h;
					N_u_prev=N_u;
					N_v_prev=N_v;

				}
				else if(param_nonlinear==2) //Linear with nonlinear advection only (valid for nondivergent nonlinear sw flows)
				{
					// First calculate the linear part
					rexiSWE.run_timestep( H, U, V, o_dt, op, simVars, param_rexi_zero_before_solving);
					//rexiSWE.run_timestep_direct_solution( H, U, V, o_dt, op, simVars );


					//Now interpolate to the the departure points
					//Departure points are set for physical space

					sampler2D.bicubic_scalar( H, posx_d, posy_d, prog_h, stag_h[0],	stag_h[1]);
					sampler2D.bicubic_scalar( U, posx_d, posy_d, prog_u, stag_u[0], stag_u[1]);
					sampler2D.bicubic_scalar( V, posx_d, posy_d, prog_v, stag_v[0], stag_v[1]);

				}

			}
			else // linear solver
			{
				rexiSWE.run_timestep( prog_h, prog_u, prog_v, o_dt,	op,	simVars, param_rexi_zero_before_solving	);
			}

			//Add forcing
			prog_h=prog_h+force_h*o_dt;
			prog_u=prog_u+force_u*o_dt;
			prog_v=prog_v+force_v*o_dt;
		}
		else if (param_timestepping_mode == 2) //Direct solution
		{
//			assert(i_max_simulation_time < 0);
			if (param_use_staggering)
			{
				std::cerr << "Direct solution on staggered grid not supported!" << std::endl;
				exit(1);
			}

			if (param_nonlinear>0)
			{
				std::cerr << "Direct solution not available for nonlinear case!" << std::endl;
				exit(1);
			}

			// Analytical solution
			assert(simVars.sim.CFL < 0);
			o_dt = -simVars.sim.CFL;
			rexiSWE.run_timestep_direct_solution(
					prog_h, prog_u, prog_v,
					-simVars.sim.CFL,
					op,
					simVars
			);
		}
		else if (param_timestepping_mode == 3)
		{   //  Implicit time step - needs check...
			assert(simVars.sim.CFL < 0);

			o_dt = -simVars.sim.CFL;
			rexiSWE.run_timestep_implicit_ts(
					prog_h, prog_u, prog_v,
					o_dt,
					op,
					simVars
			);
		}
		else if (param_timestepping_mode == 4)
		{ // Semi-Lagrangian with FD in space
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
							stag_displacement
					);

			prog_u_prev = prog_u;
			prog_v_prev = prog_v;
			prog_h_prev = prog_h;

			//Departure points are set for physical space
			//    interpolate to h grid
			sampler2D.bicubic_scalar(
					prog_h_prev,
					posx_d,
					posy_d,
					prog_h,
					stag_h[0],
					stag_h[1]
			);

		}
		else
		{
			std::cerr << "Invalid time stepping method" << std::endl;
			exit(1);
		}



		// provide information to parameters
		simVars.timecontrol.current_timestep_size = o_dt;
		simVars.timecontrol.current_simulation_time += o_dt;
		simVars.timecontrol.current_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif
	}



public:
	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{

		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();

			// Print header
			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				//if ((simVars.setup.scenario >= 0 && simVars.setup.scenario <= 4) || simVars.setup.scenario == 13)
				o_ostream << "\tDIFF_H0\tDIFF_U0\tDIFF_V0";

				if (param_compute_error && param_nonlinear==0){
					o_ostream << "\tANAL_DIFF_RMS_P\tANAL_DIFF_RMS_U\tANAL_DIFF_RMS_V";
					o_ostream << "\tANAL_DIFF_MAX_P\tANAL_DIFF_MAX_U\tANAL_DIFF_MAX_V";
				}
				o_ostream << std::endl;
			}

			//Dump  data in csv, if requested
			if (simVars.misc.output_file_name_prefix.size() > 0)
			{
				// output each time step
				if (simVars.misc.output_each_sim_seconds < 0)
					simVars.misc.output_next_sim_seconds = 0;

				if (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				{
					if (simVars.misc.output_each_sim_seconds > 0)
						while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
							simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;

					double secs = simVars.timecontrol.current_simulation_time;
					double msecs = 1000000.*(simVars.timecontrol.current_simulation_time - floor(simVars.timecontrol.current_simulation_time));
					char t_buf[256];
					sprintf(	t_buf,
								"%08d.%06d",
								(int)secs, (int)msecs
						);

					std::string ss = simVars.misc.output_file_name_prefix+"_t"+t_buf;

					prog_h.file_saveData_ascii((ss+"_h.csv").c_str());
					prog_u.file_saveData_ascii((ss+"_u.csv").c_str());
					prog_v.file_saveData_ascii((ss+"_v.csv").c_str());

					(op.diff_c_x(prog_v) - op.diff_c_y(prog_u)).file_saveData_ascii((ss+"_q.csv").c_str());
				}
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
			benchmark_diff_h = (prog_h-t0_prog_h).reduce_maxAbs() ;
			o_ostream << "\t" << benchmark_diff_h;

			// Velocity u
			benchmark_diff_u = (prog_u-t0_prog_u).reduce_maxAbs();
			o_ostream << "\t" << benchmark_diff_u;

			// Velocity v
			benchmark_diff_v = (prog_v-t0_prog_v).reduce_maxAbs();
			o_ostream << "\t" << benchmark_diff_v;
			//}

			if (param_compute_error && param_nonlinear==0)
			{
				compute_errors();

				o_ostream << "\t" << benchmark_analytical_error_rms_h;
				o_ostream << "\t" << benchmark_analytical_error_rms_u;
				o_ostream << "\t" << benchmark_analytical_error_rms_v;

				o_ostream << "\t" << benchmark_analytical_error_maxabs_h;
				o_ostream << "\t" << benchmark_analytical_error_maxabs_u;
				o_ostream << "\t" << benchmark_analytical_error_maxabs_v;
			}

			o_ostream << std::endl;
		}
	}



public:
	void compute_errors()
	{
		// Compute exact solution for linear part and compare with numerical solution

		// Initial conditions (may be in a stag grid)
		DataArray<2> t0_h = t0_prog_h;
		DataArray<2> t0_u = t0_prog_u;
		DataArray<2> t0_v = t0_prog_v;

		//Variables on unstaggered A-grid
		DataArray<2> t_h(simVars.disc.res);
		DataArray<2> t_u(simVars.disc.res);
		DataArray<2> t_v(simVars.disc.res);

		//Analytical solution at specific time on orginal grid (stag or not)
		DataArray<2> ts_h(simVars.disc.res);
		DataArray<2> ts_u(simVars.disc.res);
		DataArray<2> ts_v(simVars.disc.res);

		//The direct spectral solution can only be calculated for A grid
		t_h=t0_h;
		t_u=t0_u;
		t_v=t0_v;
		if (param_use_staggering) // Remap in case of C-grid
		{
			//remap initial condition to A grid
			//sampler2D.bilinear_scalar(t0_u, posx_a, posy_a, t_u, stag_u[0], stag_u[1]);
			//t_u=op.avg_f_x(t0_u); //equiv to bilinear
			sampler2D.bicubic_scalar(t0_u, posx_a, posy_a, t_u, stag_u[0], stag_u[1]);

			//sampler2D.bilinear_scalar(t0_v, posx_a, posy_a, t_v, stag_v[0], stag_v[1]);
			//t_v=op.avg_f_y(t0_v); //equiv to bilinear
			sampler2D.bicubic_scalar(t0_v, posx_a, posy_a, t_v, stag_v[0], stag_v[1]);
		}

		//Run exact solution for linear case
		rexiSWE.run_timestep_direct_solution(t_h, t_u, t_v, simVars.timecontrol.current_simulation_time, op, simVars);

		// Recover data in C grid using interpolations
		ts_h=t_h;
		ts_u=t_u;
		ts_v=t_v;
		if (param_use_staggering)
		{
			// Remap A grid to C grid

			//Temporary displacement for U points
			//sampler2D.bilinear_scalar(t_u, pos_x, tmp, ts_u, stag_h[0], stag_h[1]);
			//t_u=op.avg_b_x(t0_u); //equiv to bilinear
			sampler2D.bicubic_scalar(t_u, pos_x, pos_y, ts_u, stag_h[0], stag_h[1]+0.5);

			//Temporary displacement for V points
			//sampler2D.bilinear_scalar(t_v, tmp, pos_y, ts_v, stag_h[0], stag_h[1]);
			//t_v=op.avg_b_y(t0_v); //equiv to bilinear
			sampler2D.bicubic_scalar(t_v, pos_x, pos_y, ts_v, stag_h[0]+0.5, stag_h[1]);
		}

		benchmark_analytical_error_rms_h = (ts_h-prog_h).reduce_rms_quad();
		benchmark_analytical_error_rms_u = (ts_u-prog_u).reduce_rms_quad();
		benchmark_analytical_error_rms_v = (ts_v-prog_v).reduce_rms_quad();

		benchmark_analytical_error_maxabs_h = (ts_h-prog_h).reduce_maxAbs();
		benchmark_analytical_error_maxabs_u = (ts_u-prog_u).reduce_maxAbs();
		benchmark_analytical_error_maxabs_v = (ts_v-prog_v).reduce_maxAbs();

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
		const DataArray<2>* data;
		const char *description;
	};

	/**
	 * Arrays for online visualisation and their textual description
	 */
	VisStuff vis_arrays[4] =
	{
			{&prog_h,	"h"},
			{&prog_u,	"u"},
			{&prog_v,	"v"},
			{&eta,		"eta"}
	};



	void vis_get_vis_data_array(
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		if (simVars.misc.vis_id < 0)
		{
			DataArray<2> t_h = t0_prog_h;
			DataArray<2> t_u = t0_prog_u;
			DataArray<2> t_v = t0_prog_v;

			rexiSWE.run_timestep_direct_solution(
					t_h, t_u, t_v,
					simVars.timecontrol.current_simulation_time,
					op,
					simVars
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

			case -4:					// beta plane
				vis = beta_plane;
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
			prog_h.file_saveData_ascii("swe_rexi_dump_h.csv");
			prog_u.file_saveData_ascii("swe_rexi_dump_u.csv");
			prog_v.file_saveData_ascii("swe_rexi_dump_v.csv");
			break;

		case 'C':
			// dump data arrays to VTK
			prog_h.file_saveData_vtk("swe_rexi_dump_h.vtk", "Height");
			prog_u.file_saveData_vtk("swe_rexi_dump_u.vtk", "U-Velocity");
			prog_v.file_saveData_vtk("swe_rexi_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_h.file_loadData("swe_rexi_dump_h.csv", simVars.setup.input_data_binary);
			prog_u.file_loadData("swe_rexi_dump_u.csv", simVars.setup.input_data_binary);
			prog_v.file_loadData("swe_rexi_dump_v.csv", simVars.setup.input_data_binary);
			break;
		}
	}


	bool instability_detected()
	{
		return !(	prog_h.reduce_all_finite() &&
					prog_u.reduce_all_finite() &&
					prog_v.reduce_all_finite()
				);
	}


#if SWEET_PARAREAL

	/******************************************************
	 ******************************************************
	 *       ************** PARAREAL **************
	 ******************************************************
	 ******************************************************/

	DataArray<2> _parareal_data_start_h, _parareal_data_start_u, _parareal_data_start_v;
	Parareal_Data_DataArrays<3> parareal_data_start;

	DataArray<2> _parareal_data_fine_h, _parareal_data_fine_u, _parareal_data_fine_v;
	Parareal_Data_DataArrays<3> parareal_data_fine;

	DataArray<2> _parareal_data_coarse_h, _parareal_data_coarse_u, _parareal_data_coarse_v;
	Parareal_Data_DataArrays<3> parareal_data_coarse;

	DataArray<2> _parareal_data_output_h, _parareal_data_output_u, _parareal_data_output_v;
	Parareal_Data_DataArrays<3> parareal_data_output;

	DataArray<2> _parareal_data_error_h, _parareal_data_error_u, _parareal_data_error_v;
	Parareal_Data_DataArrays<3> parareal_data_error;

	double timeframe_start = -1;
	double timeframe_end = -1;

	bool output_data_valid = false;

	void parareal_setup()
	{
		{
			DataArray<2>* data_array[3] = {&_parareal_data_start_h, &_parareal_data_start_u, &_parareal_data_start_v};
			parareal_data_start.setup(data_array);
		}

		{
			DataArray<2>* data_array[3] = {&_parareal_data_fine_h, &_parareal_data_fine_u, &_parareal_data_fine_v};
			parareal_data_fine.setup(data_array);
		}

		{
			DataArray<2>* data_array[3] = {&_parareal_data_coarse_h, &_parareal_data_coarse_u, &_parareal_data_coarse_v};
			parareal_data_coarse.setup(data_array);
		}

		{
			DataArray<2>* data_array[3] = {&_parareal_data_output_h, &_parareal_data_output_u, &_parareal_data_output_v};
			parareal_data_output.setup(data_array);
		}

		{
			DataArray<2>* data_array[3] = {&_parareal_data_error_h, &_parareal_data_error_u, &_parareal_data_error_v};
			parareal_data_error.setup(data_array);
		}

		// use REXI
		rexiSWE.setup(
				param_rexi_h,
				param_rexi_m,
				param_rexi_l,

				simVars.disc.res,
				simVars.sim.domain_size,
				param_rexi_half,
				param_rexi_use_spectral_differences_for_complex_array,
				param_rexi_helmholtz_solver_id,
				param_rexi_helmholtz_solver_eps
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

		// cast to pararealDataArray stuff
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

		// run implicit time step
//		assert(i_max_simulation_time < 0);
//		assert(simVars.sim.CFL < 0);

		rexiSWE.run_timestep_implicit_ts(
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
		if (simVars.parareal.verbosity > 2)
			std::cout << "compute_output_data()" << std::endl;

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

		if (param_compute_error && param_nonlinear==0){
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
		Parareal_Data_DataArrays<3>& data = (Parareal_Data_DataArrays<3>&)i_data;

		std::ostringstream ss;
		ss << "output_iter" << iteration_id << "_slice" << time_slice_id << ".vtk";

		std::string filename = ss.str();

//		std::cout << "filename: " << filename << std::endl;
		data.data_arrays[0]->file_saveData_vtk(filename.c_str(), filename.c_str());
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
		{
				std::cerr << "MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye..." << std::endl;
				exit(-1);
		}
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

	NUMABlockAlloc::setup();

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"rexi-h",			///Rexi parameters
			"rexi-m",
			"rexi-l",
			"rexi-half",
			"timestepping-mode",			///Method to be used
			"compute-error",
			"staggering",						/// FD A-C grid
			"use-specdiff-for-complex-array",	/// use finite differences for complex array
			"rexi-helmholtz-solver-id",		/// use iterative solver for REXI
			"rexi-helmholtz-solver-eps",		/// error threshold for solver
			"boundary-id",
			"rexi-zero-before-solving",
			"nonlinear",    /// form of equations
			"initial-freq-x-mul",		// frequency multipliers for special scenario setup
			"initial-freq-y-mul",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 0.2;  // REXI h
	simVars.bogus.var[1] = 256;	// M
	simVars.bogus.var[2] = 0;	// L = 0: default - this gives L=11
	simVars.bogus.var[3] = 1;	// param_rexi_half
	simVars.bogus.var[4] = 0;	// timestepping mode - default is RK
	simVars.bogus.var[5] = 0;	// compute error - default no
	simVars.bogus.var[6] = 0;	// stag - default A grid
	simVars.bogus.var[7] = 1;	// Use spec diff per default for complex array
	simVars.bogus.var[8] = 0;	// Use spectral solver id
	simVars.bogus.var[9] = 1e-7;	// Error threshold
	simVars.bogus.var[10] = 0; 	//boundary
	simVars.bogus.var[11] = 1;	// zero rexi
	simVars.bogus.var[12] = 0;	// nonlinear
	simVars.bogus.var[13] = 0;
	simVars.bogus.var[14] = 0;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--rexi-h [h-value]	h-sampling distance for REXI, default:0.2" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--rexi-m [m-value]	M-value for REXI (related to number of poles), default:256" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--rexi-half [0/1]	Reduce rexi computations to its half, default:1 " << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--timestepping-mode [0/1/2/3/4]	Timestepping method to use" << std::endl;
		std::cout << "	                            0: RKn with Finite-difference (default)" << std::endl;
		std::cout << "	                            1: REXI" << std::endl;
		std::cout << "	                            2: Direct solution in spectral space" << std::endl;
		std::cout << "	                            3: Implicit Finite-difference (needs checking)" << std::endl;
		std::cout << "	                            4: Semi-Lagrangian with Finite-difference" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--compute-error [0/1]	Compute the errors" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--staggering [0/1]		Use staggered grid" << std::endl;
		std::cout << std::endl;
		std::cout << "	--use-specdiff-for-complex-array [0/1]	Controls the discretization of the derivative operations for Complex2DArrays:" << std::endl;
		std::cout << "                                      0: Finite-difference derivatives" << std::endl;
		std::cout << "                                      1: Spectral derivatives (default)" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-helmholtz-solver-id [int]	Use iterative solver for REXI" << std::endl;
		std::cout << "	--rexi-helmholtz-solver-eps [err]	Error threshold for iterative solver" << std::endl;
		std::cout << std::endl;
		std::cout << "	--boundary-id [0,1,...]	    Boundary id" << std::endl;
		std::cout << "                              0: no boundary (default)" << std::endl;
		std::cout << "                              1: centered box" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-zero-before-solving [0/1]	Zero the solution for the iterative solver (default=0)" << std::endl;
		std::cout << std::endl;
		std::cout << "	--nonlinear [0/1/2]	   Form of equations:" << std::endl;
		std::cout << "						     0: Linear SWE (default)" << std::endl;
		std::cout << "						     1: Full nonlinear SWE" << std::endl;
		std::cout << "						     2: Linear SWE + nonlinear advection only (needs -H to be set)" << std::endl;
		std::cout << std::endl;


#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}

	//Rexi parameters
	param_rexi_h = simVars.bogus.var[0];
    param_rexi_m = simVars.bogus.var[1];
	param_rexi_l = simVars.bogus.var[2];
	param_rexi_half = simVars.bogus.var[3];

	// Method to be used
	param_timestepping_mode = simVars.bogus.var[4];

	//Calculate error flag
	param_compute_error = simVars.bogus.var[5];

	// C- grid flag
	param_use_staggering = simVars.bogus.var[6];

	// Linear solver
	param_rexi_use_spectral_differences_for_complex_array = simVars.bogus.var[7];
    param_rexi_helmholtz_solver_id = simVars.bogus.var[8];
	param_rexi_helmholtz_solver_eps = simVars.bogus.var[9];

	//Boundary
	param_boundary_id = simVars.bogus.var[10];

	param_rexi_zero_before_solving = simVars.bogus.var[11];

	// Linear vs nonlinear swe
	param_nonlinear = simVars.bogus.var[12];

	param_initial_freq_x_mul = simVars.bogus.var[13];
	param_initial_freq_y_mul = simVars.bogus.var[14];


	//Print header
	std::cout << std::endl;
	if(param_nonlinear == 1)
		std::cout << "Solving full nonlinear SW equations" << std::endl;
	else
	{
		if(param_nonlinear == 2)
				std::cout << "Solving linear SWE with nonlinear advection" << std::endl;
		else
			std::cout << "Solving linear SW equations" << std::endl;
	}

	std::cout << "-----------------------------" << std::endl;
	std::cout << "Method to be used (timestepping_mode): "; // << param_timestepping_mode << std::endl;
	switch(param_timestepping_mode)
	{
		case 0:
				std::cout << " 0: Finite-difference " << std::endl; break;
		case 1:
				std::cout << " 1: REXI" << std::endl; break;
		case 2:
				std::cout << " 2: Direct solution in spectral space" << std::endl; break;
		case 3:
				std::cout << " 3: Implicit method - needs checking!" << std::endl; break;
		case 4:
				std::cout << " 4: Semi-Lag with FD - needs checking!" << std::endl; break;
		default:
			std::cerr << "Timestepping unknowkn" << std::endl;
			return -1;
	}
	std::cout << "Staggered grid? " << param_use_staggering << std::endl;
	std::cout << "Computing error? " << param_compute_error << std::endl;
	std::cout << "Verbosity: " << simVars.misc.verbosity << std::endl;
	std::cout << "Parareal: " << SWEET_PARAREAL << std::endl;

	std::ostringstream buf;
	buf << std::setprecision(14);


#if SWEET_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// only start simulation and time stepping for first rank
	if (rank == 0)
#endif
	{
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

			//Main time loop
			while(true)
			{
				//Output data
				if (simVars.misc.verbosity > 1)
				{
					simulationSWE->timestep_output(buf);

					std::string output = buf.str();
					buf.str("");

					//This is an output printed on screen or buffered to files if > used
					std::cout << output;

					//This is an output only printed on screen, not buffered to files if > used
					if (simVars.misc.verbosity > 2)
						std::cerr << output;
				}

				//Stop simulation if requested
				if (simulationSWE->should_quit())
					break;

				//Main call for timestep run
				simulationSWE->run_timestep();

				//Instability
				if (simulationSWE->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}
			}

			//Stop counting time
			time.stop();

			double seconds = time();

			//End of run output results
			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;
			if (param_timestepping_mode != 0)
				std::cout << "REXI alpha.size(): " << simulationSWE->rexiSWE.rexi.alpha.size() << std::endl;

			if (simVars.misc.verbosity > 0)
			{
				std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((simVars.diag.total_energy-diagnostics_energy_start)/diagnostics_energy_start) << std::endl;
				std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
				std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((simVars.diag.total_potential_enstrophy-diagnostics_potential_entrophy_start)/diagnostics_potential_entrophy_start) << std::endl;

				if (param_compute_error)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark_diff_h << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
				}
			}

			if (param_compute_error && param_nonlinear==0)
			{
				simulationSWE->compute_errors();
				std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << simulationSWE->benchmark_analytical_error_rms_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << simulationSWE->benchmark_analytical_error_rms_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << simulationSWE->benchmark_analytical_error_rms_v << std::endl;

				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << simulationSWE->benchmark_analytical_error_maxabs_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << simulationSWE->benchmark_analytical_error_maxabs_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << simulationSWE->benchmark_analytical_error_maxabs_v << std::endl;
			}

			delete simulationSWE;
		}
	}
#if SWEET_MPI
	else
	{
		if (param_timestepping_mode == 1)
		{
			RexiSWE rexiSWE;

			/*
			 * Setup our little dog REXI
			 */
			rexiSWE.setup(
					param_rexi_h,
					param_rexi_m,
					param_rexi_l,
					simVars.disc.res,
					simVars.sim.domain_size,
					param_rexi_half,
					param_rexi_use_spectral_differences_for_complex_array,
					param_rexi_helmholtz_solver_id,
					param_rexi_helmholtz_solver_eps
				);

			bool run = true;

			DataArray<2> prog_h(simVars.disc.res);
			DataArray<2> prog_u(simVars.disc.res);
			DataArray<2> prog_v(simVars.disc.res);

			Operators2D op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

			MPI_Barrier(MPI_COMM_WORLD);

			while (run)
			{
				// REXI time stepping
				run = rexiSWE.run_timestep(
						prog_h, prog_u, prog_v,
						-simVars.sim.CFL,
						op,
						simVars,
						param_rexi_zero_before_solving
				);
			}
		}
	}
#endif


#if SWEET_MPI
	if (param_timestepping_mode > 0)
	{
		// synchronize REXI
		if (rank == 0)
			RexiSWE::MPI_quitWorkers(simVars.disc.res);
	}

	MPI_Finalize();
#endif

	return 0;
}
