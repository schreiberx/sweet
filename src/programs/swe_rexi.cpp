/*
* SWM with nonlinear part using REXI test
*
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
#include "rexiswe/RexiSWE.hpp"

#ifndef SWEET_MPI
#	define SWEET_MPI 1
#endif

#if SWEET_MPI
#	include <mpi.h>
#endif

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
bool param_nonlinear;

//Main simulation class
class SimulationSWE
{
public:
	// Prognostic variables
	DataArray<2> prog_h, prog_u, prog_v;
	DataArray<2> prog_u_prev, prog_v_prev;

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

	// Arrival points for semi-lag
	DataArray<2> posx_a, posy_a;
	DataArray<2> *pos_a[2];

	//Staggering displacement array (use 0.5 for each displacement)
	// [0] - delta x of u variable
	// [1] - delta y of u variable
	// [2] - delta x of v variable
	// [3] - delta y of v variable
	// For C grid use {-0.5,0,0,-0.5}
	// Default - A grid (zeros)
	double stag_displacement[4] = {0,0,0,0};

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
	SimulationSWE()	:
	// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_h(simVars.disc.res),
		prog_u(simVars.disc.res),
		prog_v(simVars.disc.res),

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

		posx_a(simVars.disc.res),
		posy_a(simVars.disc.res),

		//h_t(simVars.disc.res),

		// Initialises operators
		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
	{
		// Calls initialisation of the run (e.g. sets u, v, h)
		reset();
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

		sampler2D.setup(simVars.sim.domain_size, simVars.disc.res);

		//Check if input parameters are adequate for this simulation
		if (param_use_staggering && simVars.disc.use_spectral_basis_diffs)
		{
			std::cerr << "Staggering and spectral basis not supported!" << std::endl;
			exit(1);
		}

		if (param_use_staggering && param_timestepping_mode != 0)
		{
			std::cerr << "Staggering only supported for standard time stepping mode 0!" << std::endl;
			exit(1);
		}

		if (param_use_staggering)
		{
			/*
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
			 *                            udx  udy  vdx   vdy   */
			double stag_displacement[4] = {0.0, 0.5, 0.5, 0.0};
		}


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

						prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
						t0_prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));

						//PXT - Why is t0 here??? This makes the error calculation wrong for the FD case with C grid
						//t0_prog_u.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
						//t0_prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
					}

					{
						// u space
						double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						prog_u.set(j,i, SWEValidationBenchmarks::return_u(simVars, x, y));
						t0_prog_u.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
					}

					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
						t0_prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
					}
				}
				else // A-Grid (colocated grid)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

					prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
					prog_u.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
					prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));

					t0_prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
					t0_prog_u.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
					t0_prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
				}
			}
		}


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
		if (param_timestepping_mode == 1 || param_timestepping_mode == 3)
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

		//if(param_use_staggering){
		//	if (simVars.sim.beta != 0)
		//	{
		//		q = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + beta_plane) / prog_h;
		//	}
		//	else
		//	{
		//		q = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + simVars.sim.f0) / prog_h;
		//	}
		//}


		// mass
		simVars.diag.total_mass = prog_h.reduce_sum_quad() * normalization;

		// energy
		simVars.diag.total_energy = 0.5*((
				prog_h*prog_h +
				prog_h*prog_u*prog_u +
				prog_h*prog_v*prog_v
			).reduce_sum_quad()) * normalization;

		// potential enstropy
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
				double limit_speed = std::min(simVars.disc.cell_size[0]/i_u.reduce_maxAbs(), simVars.disc.cell_size[1]/i_v.reduce_maxAbs());
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
				double limit_gh = std::min(simVars.disc.cell_size[0], simVars.disc.cell_size[1])/std::sqrt(simVars.sim.g*i_h.reduce_maxAbs());

				if (simVars.misc.verbosity > 2)
					std::cerr << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

				o_dt = simVars.sim.CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);
			}
		}

		// A- grid method
		if (!param_use_staggering)
		{
			if(param_nonlinear){ //nonlinear case
				std::cout << "Only linear swe are setup for unstaggered grids " << std::endl;
				exit(1);
			}

			boundary_action();

			/*
			 * linearized non-conservative (advective) formulation:
			 *
			 * h_t = -h0*u_x - h0*v_y
			 * u_t = -g * h_x + f*v
			 * v_t = -g * h_y - f*u
			 */
			o_u_t = -simVars.sim.g*op.diff_c_x(i_h) + simVars.sim.f0*i_v;
			o_v_t = -simVars.sim.g*op.diff_c_y(i_h) - simVars.sim.f0*i_u;


			if (simVars.sim.viscosity != 0)
			{
				o_u_t -= op.diffN_x(i_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
				o_v_t -= op.diffN_y(i_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			}

			boundary_action();

			if (!simVars.disc.timestepping_leapfrog_like_update)
			{
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
			}
			else
			{
				/*
				 * a kind of leapfrog:
				 *
				 * We use the hew v and u values to compute the update for p
				 *
				 * compute updated u and v values without using it
				 */
				if (!simVars.disc.timestepping_up_and_downwinding)
				{
					// recompute U and V

					// update based on new u and v values
					o_h_t = -op.diff_c_x(
								simVars.setup.h0*(i_u+o_dt*o_u_t)
							)
							-op.diff_c_y(
									simVars.setup.h0*(i_v+o_dt*o_v_t)
							);
				}
				else
				{
					// update based on new u and v values
					compute_upwinding_P_updates(
							i_h,
							i_u+o_dt*o_u_t,
							i_v+o_dt*o_v_t,
							o_h_t
						);
				}
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
			if(param_nonlinear) //nonlinear case
			{
				U = op.avg_b_x(i_h)*i_u;
				V = op.avg_b_y(i_h)*i_v;
			}
			else // linear case
			{
				U = simVars.setup.h0*i_u;
				V = simVars.setup.h0*i_v;
			}

			if(param_nonlinear) //nonlinear case
				H = simVars.sim.g*i_h + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));
			else //linear case
				H = simVars.sim.g*i_h;// + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));

			if(param_nonlinear) //nonlinear case
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
			if (!simVars.disc.timestepping_leapfrog_like_update)
			{
				if (!simVars.disc.timestepping_up_and_downwinding)
				{
					// standard update
					o_h_t = -op.diff_f_x(U) - op.diff_f_y(V);
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
			else
			{
				/*
				 * a kind of leapfrog (this is not conserving anything! don't use it for production runs!):
				 *
				 * We use the new v and u values to compute the update for p
				 *
				 * compute updated u and v values without using it
				 */
				if (!simVars.disc.timestepping_up_and_downwinding)
				{
					// recompute U and V
					U = op.avg_b_x(i_h)*(i_u+o_dt*o_u_t);
					V = op.avg_b_y(i_h)*(i_v+o_dt*o_v_t);

					// update based on new u and v values
					o_h_t = -op.diff_f_x(U) - op.diff_f_y(V);
				}
				else
				{
					// update based on new u and v values
					compute_upwinding_P_updates(
							i_h,
							i_u+o_dt*o_u_t,
							i_v+o_dt*o_v_t,
							o_h_t
						);
				}
			}

#if 0
			if (simVars.sim.potential_viscosity != 0)
				o_h_t -= op.diff2(i_h)*simVars.sim.potential_viscosity;

			if (simVars.sim.potential_hyper_viscosity != 0)
				o_h_t -= op.diff4(i_h)*simVars.sim.potential_hyper_viscosity;
#endif
		}
	}


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

	void run_timestep()
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
					&SimulationSWE::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					prog_h, prog_u, prog_v,
					o_dt,
					simVars.timecontrol.current_timestep_size,
					simVars.disc.timestepping_runge_kutta_order,
					simVars.timecontrol.current_simulation_time,
					simVars.timecontrol.max_simulation_time
				);
		}
		else if (param_timestepping_mode == 1)
		{
			// REXI time stepping for linear eq
			if (param_nonlinear)
			{
				o_dt = -simVars.sim.CFL;
			}
			else // linear solver
			{
				o_dt = -simVars.sim.CFL;
				rexiSWE.run_timestep(
						prog_h, prog_u, prog_v,
						o_dt,
						op,
						simVars,
						param_rexi_zero_before_solving
				);
			}
		}
		else if (param_timestepping_mode == 2)
		{
			if (param_use_staggering)
			{
				std::cerr << "Direct solution on staggered grid not supported!" << std::endl;
				exit(1);
			}

			if (param_nonlinear)
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
		{
			assert(simVars.sim.CFL < 0);
			// Analytical solution
			o_dt = -simVars.sim.CFL;
			rexiSWE.run_timestep_implicit_ts(
					prog_h, prog_u, prog_v,
					o_dt,
					op,
					simVars
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

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
					o_ostream << "\tDIFF_P0\tDIFF_U0\tDIFF_V0";

				if (param_compute_error){
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
			o_ostream << std::setprecision(6) << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;

			// Print max abs difference of vars to initial conditions (this gives the error in steady state cases)
			if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
			{
				// Height
				benchmark_diff_h = (prog_h-t0_prog_h).reduce_maxAbs() ;
				o_ostream << "\t" << benchmark_diff_h;

				// Velocity u
				benchmark_diff_u = (prog_u-t0_prog_u).reduce_maxAbs();
				o_ostream << "\t" << benchmark_diff_u;

				// Velocity v
				benchmark_diff_v = (prog_v-t0_prog_v).reduce_maxAbs();
				o_ostream << "\t" << benchmark_diff_v;
			}

			if (param_compute_error && !param_nonlinear)
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
		// Initial conditions (may be in a stag grid)
		DataArray<2> t0_h = t0_prog_h;
		DataArray<2> t0_u = t0_prog_u;
		DataArray<2> t0_v = t0_prog_v;
		//Variables for remapping t0/from staggered grid
		DataArray<2> t_h(simVars.disc.res);
		DataArray<2> t_u(simVars.disc.res);
		DataArray<2> t_v(simVars.disc.res);

		if(param_nonlinear)
		{
			std::cout << "Warning: Exact solution not possible in general for nonlinear swe. Using exact solution for linear case instead." << std::endl;
		}

		//std::cout << std::endl;
		//std::cout << t_v << std::endl;

		//The direct spectral solution can only be calculated for A grid
		if (param_use_staggering)
		{
			std::cout << "Warning: Staggered data being interpolated to A-grid for exact linear solution" << std::endl;
			//remap initial condition to A grid
			sampler2D.remap_gridC2A(t0_u, t0_v, t_u, t_v, 0);
		}
		//std::cout << std::endl;
		//std::cout << t_h << std::endl;

		//Run exact solution for linear case
		rexiSWE.run_timestep_direct_solution(t_h, t_u, t_v, simVars.timecontrol.current_simulation_time, op, simVars);


		if (param_use_staggering)
		{
			//remap exact solution to C grid - TODO
		}

		benchmark_analytical_error_rms_h = (t_h-prog_h).reduce_rms_quad();
		if (!param_use_staggering)
		{
			benchmark_analytical_error_rms_u = (t_u-prog_u).reduce_rms_quad();
			benchmark_analytical_error_rms_v = (t_v-prog_v).reduce_rms_quad();
		}

		benchmark_analytical_error_maxabs_h = (t_h-prog_h).reduce_maxAbs();
		if (!param_use_staggering)
		{
			benchmark_analytical_error_maxabs_u = (t_u-prog_u).reduce_maxAbs();
			benchmark_analytical_error_maxabs_v = (t_v-prog_v).reduce_maxAbs();
		}
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
			//PXT - Why do this and not simply use t0_prog?
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
				vis = t_h-prog_h;	// difference
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
};



int main(int i_argc, char *i_argv[])
{

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
			"rexi-h",
			"rexi-m",
			"rexi-l",
			"rexi-half",
			"timestepping-mode",
			"compute-error",
			"staggering",
			"use-specdiff-for-complex-array",	/// use finite differences for complex array
			"rexi-helmholtz-solver-id",		/// use iterative solver for REXI
			"rexi-helmholtz-solver-eps",		/// error threshold for solver
			"boundary-id",
			"rexi-zero-before-solving",
			"nonlinear",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 0.2;
	simVars.bogus.var[1] = 256;	// M
	simVars.bogus.var[2] = 0;	// L = 0: default - this gives L=11
	simVars.bogus.var[3] = 1;	// param_rexi_half
	simVars.bogus.var[4] = 0;	// timestepping mode
	simVars.bogus.var[5] = 0;	// compute error
	simVars.bogus.var[6] = 0;	// stag
	simVars.bogus.var[7] = 1;	// Use spec diff per default for complex array
	simVars.bogus.var[8] = 0;	// Use spectral solver id
	simVars.bogus.var[9] = 1e-7;	// Error threshold
	simVars.bogus.var[10] = 0; 	//boundary
	simVars.bogus.var[11] = 1;	// zero rexi
	simVars.bogus.var[12] = 0;	// nonlinear

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
		std::cout << "	--timestepping-mode [0/1/2]	Timestepping method to use" << std::endl;
		std::cout << "	                            0: RKn with Finite-difference (default)" << std::endl;
		std::cout << "	                            1: REXI" << std::endl;
		std::cout << "	                            2: Direct solution in spectral space" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--compute-error [0/1]	Compute the errors" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--staggering [0/1]		Use staggered grid" << std::endl;
		std::cout << std::endl;
		std::cout << "	--use-specdiff-for-complex-array=[0/1]	Controls the discretization of the derivative operations for Complex2DArrays:" << std::endl;
		std::cout << "                                      0: Finite-difference derivatives" << std::endl;
		std::cout << "                                      1: Spectral derivatives (default)" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-helmholtz-solver-id=[int]	Use iterative solver for REXI" << std::endl;
		std::cout << "	--rexi-helmholtz-solver-eps=[err]	Error threshold for iterative solver" << std::endl;
		std::cout << std::endl;
		std::cout << "	--boundary-id=[0,1,...]	    Boundary id" << std::endl;
		std::cout << "                              0: no boundary" << std::endl;
		std::cout << "                              1: centered box" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-zero-before-solving=[0/1]	Zero the solution for the iterative solver" << std::endl;
		std::cout << std::endl;
		std::cout << "	--nonlinear=[0/1]	   0-use linear equations (default), 1-use nonlinear eq" << std::endl;
		std::cout << std::endl;
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

	//Print header
	std::cout << std::endl;
	if(param_nonlinear)
		std::cout << "Solving nonlinear SW equations" << std::endl;
	else
		std::cout << "Solving linear SW equations" << std::endl;

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
		default:
			std::cerr << "Timestepping unknowkn" << std::endl;
			return -1;
	}
	std::cout << "Staggered grid? " << param_use_staggering << std::endl;
	std::cout << "Computing error? " << param_compute_error << std::endl;
	std::cout << "Verbosity: " << simVars.misc.verbosity << std::endl;

	SimulationSWE *simulationSWE = new SimulationSWE;

	std::ostringstream buf;
	buf << std::setprecision(14);

#if SWEET_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// only start simulation and time stepping for first rank
	if (rank == 0)
#endif
	{
#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			VisSweet<SimulationSWE> visSweet(simulationSWE);
		}
		else
#endif
		{
			//Setting initial conditions and workspace - in case there is no GUI
			simulationSWE->reset();

			//Time counter
			Stopwatch time;

			//Diagnostic measures at initial stage
			double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

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

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark_diff_h << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
				}
			}

			if (param_compute_error)
			{
				simulationSWE->compute_errors();

				std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << simulationSWE->benchmark_analytical_error_rms_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << simulationSWE->benchmark_analytical_error_rms_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << simulationSWE->benchmark_analytical_error_rms_v << std::endl;

				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << simulationSWE->benchmark_analytical_error_maxabs_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << simulationSWE->benchmark_analytical_error_maxabs_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << simulationSWE->benchmark_analytical_error_maxabs_v << std::endl;
			}
		}

		delete simulationSWE;
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
