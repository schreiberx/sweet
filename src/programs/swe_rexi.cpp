/*
* SWM with nonlinear part using REXI test
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
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include "rexi/RexiSWE.hpp"

#ifndef SWEET_MPI
#	define SWEET_MPI 1
#endif

#if SWEET_MPI
#	include <mpi.h>
#endif

SimulationVariables simVars;

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

double param_initial_freq_x_mul;
double param_initial_freq_y_mul;

int param_boundary_id;
bool param_nonlinear;

class SimulationSWE
{
public:
	DataArray<2> prog_h, prog_u, prog_v;
	DataArray<2> eta, q;
	DataArray<2> tmp;

	// temporary variables;
	DataArray<2> tmp0, tmp1, tmp2, tmp3;

	DataArray<2> boundary_mask;
	DataArray<2> boundary_mask_inv;

	// initial values for comparison with analytical solution
	DataArray<2> t0_prog_h, t0_prog_u, t0_prog_v;

	Operators2D op;

	TimesteppingRK timestepping;

	int last_timestep_nr_update_diagnostics = -1;

	double benchmark_diff_h;
	double benchmark_diff_u;
	double benchmark_diff_v;

	double benchmark_analytical_error_rms_h;
	double benchmark_analytical_error_rms_u;
	double benchmark_analytical_error_rms_v;

	double benchmark_analytical_error_maxabs_h;
	double benchmark_analytical_error_maxabs_u;
	double benchmark_analytical_error_maxabs_v;

	RexiSWE rexiSWE;


public:
	SimulationSWE()	:
		prog_h(simVars.disc.res),
		prog_u(simVars.disc.res),
		prog_v(simVars.disc.res),

		eta(simVars.disc.res),
		q(simVars.disc.res),
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

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
	{
		reset();
	}


	void reset()
	{
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

		prog_h.set_all(simVars.setup.h0);
		prog_u.set_all(0);
		prog_v.set_all(0);
		boundary_mask.set_all(0);

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

#if 0
		if (simVars.setup.scenario == -1)
		{
			std::cout << "Setting up steady state" << std::endl;

			if (simVars.sim.f0 == 0)
			{
				std::cerr << "Coriolis frequency f is set to 0" << std::endl;
				exit(1);
			}

			if (std::isinf(param_rexi_h))
				std::cerr << "ERROR: Steady state scenario here only makes sense for linear solver only" << std::endl;

			/**
			 * setup steady state
			 */
			if (param_use_staggering)
			{
				std::cerr << "NOT SUPPORTED YET" << std::endl;
				exit(1);
			}
			else
			{
				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
				{
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						// Gaussian
						double dx = x-simVars.setup.coord_x*simVars.sim.domain_size[0];
						double dy = y-simVars.setup.coord_y*simVars.sim.domain_size[1];

						double radius = simVars.setup.radius_scale*sqrt((double)simVars.sim.domain_size[0]*(double)simVars.sim.domain_size[0]+(double)simVars.sim.domain_size[1]*(double)simVars.sim.domain_size[1]);
						dx /= radius;
						dy /= radius;

#if 1
						double sx = x/simVars.sim.domain_size[0];
						double sy = y/simVars.sim.domain_size[1];

						double foo = std::sin(sx*2.0*M_PIl)*std::sin(sy*2.0*M_PIl);

						double barx = std::cos(sx*2.0*M_PIl)*std::sin(sy*2.0*M_PIl)*2.0*M_PIl/simVars.sim.domain_size[0];
						double bary = std::sin(sx*2.0*M_PIl)*std::cos(sy*2.0*M_PIl)*2.0*M_PIl/simVars.sim.domain_size[1];

						double h = simVars.setup.h0+foo;
						double u = -simVars.sim.g/simVars.sim.f0*bary;
						double v = simVars.sim.g/simVars.sim.f0*barx;
#else
						double foo = std::exp(-50.0*(dx*dx + dy*dy));
						double h = simVars.setup.h0+foo;
						double u = -simVars.sim.g/simVars.sim.f0*(-50.0*2.0*dy)*foo;
						double v = simVars.sim.g/simVars.sim.f0*(-50.0*2.0*dx)*foo;
#endif

						prog_h.set(j, i, h);
						prog_u.set(j, i, u);
						prog_v.set(j, i, v);

						t0_prog_h.set(j, i, h);
						t0_prog_u.set(j, i, u);
						t0_prog_v.set(j, i, v);
					}
				}
			}
		}
		else
#endif
		if (param_initial_freq_x_mul != -1)
		{
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					double h_dx, h_dy;
					double u_dx, u_dy;
					double v_dx, v_dy;

					if (param_use_staggering)
					{
						{
							// h
							double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

							h_dx = x/simVars.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
							h_dy = y/simVars.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
						}

						{
							// u space
							double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

							u_dx = x/simVars.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
							u_dy = y/simVars.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
						}

						{
							// v space
							double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

							v_dx = x/simVars.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
							v_dy = y/simVars.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;
						}
					}
					else
					{
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						h_dx = x/simVars.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
						h_dy = y/simVars.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;

						u_dx = h_dx;
						u_dy = h_dy;

						v_dx = h_dx;
						v_dy = h_dy;
					}

					//
					// note, that cos*cos and cos*sin is a function which cannot be directly represented in spectral space
					//
					double h = std::sin(2.0*h_dx)*std::cos(2.0*h_dy) - (1.0/5.0)*std::cos(2.0*h_dx)*std::sin(4.0*h_dy) + simVars.setup.h0;
					double u = std::cos(4.0*u_dx)*std::cos(2.0*u_dy);
					double v = std::cos(2.0*v_dx)*std::cos(4.0*v_dy);

					prog_h.set(j, i, h);
					prog_u.set(j, i, u);
					prog_v.set(j, i, v);

					// REXI initial conditions given by values on non-staggered grid
					double t0_h = std::sin(2.0*h_dx)*std::cos(2.0*h_dy) - (1.0/5.0)*std::cos(2.0*h_dx)*std::sin(4.0*h_dy) + simVars.setup.h0;
					double t0_u = std::cos(4.0*h_dx)*std::cos(2.0*h_dy);
					double t0_v = std::cos(2.0*h_dx)*std::cos(4.0*h_dy);
					t0_prog_h.set(j, i, t0_h);
					t0_prog_u.set(j, i, t0_u);
					t0_prog_v.set(j, i, t0_v);
				}
			}
		}
		else
		{
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					if (param_use_staggering)
					{
						{
							// h
							double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

							prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));

							t0_prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
							t0_prog_u.set(j,i, SWEValidationBenchmarks::return_u(simVars, x, y));
							t0_prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
						}

						{
							// u space
							double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

							prog_u.set(j,i, SWEValidationBenchmarks::return_u(simVars, x, y));
						}

						{
							// v space
							double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

							prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
						}
					}
					else
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
		}


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

		if (simVars.misc.gui_enabled)
			timestep_output();
	}



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
		simVars.diag.total_energy = 0.5*(
				prog_h*prog_h +
				prog_h*prog_u*prog_u +
				prog_h*prog_v*prog_v
			).reduce_sum_quad() * normalization;

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

		// same a above, but formulated in a finite-difference style
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
//		assert(param_boundary_id != 0);

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
		if (!param_use_staggering)
		{
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
					double hx = simVars.disc.cell_size[0];
					double hy = simVars.disc.cell_size[1];

					double limit_speed = std::min(hx/i_u.reduce_maxAbs(), hy/i_v.reduce_maxAbs());

					// limit by re
					double limit_visc = std::numeric_limits<double>::infinity();

					// limit by gravitational acceleration
					double limit_gh = std::min(simVars.disc.cell_size[0], simVars.disc.cell_size[1])/std::sqrt(simVars.sim.g*i_h.reduce_maxAbs());

					if (simVars.misc.verbosity > 2)
						std::cerr << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

					o_dt = simVars.sim.CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);
				}
			}

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
		else
		{
			boundary_action();

			// STAGGERED GRID

			DataArray<2>& U = tmp0;
			DataArray<2>& V = tmp1;
			DataArray<2>& q = tmp2;
			DataArray<2>& H = tmp3;

			/*
			 * Note, that this grid does not follow the formulation
			 * in the paper of Robert Sadourny, but looks as follows:
			 *
			 *             ^
			 *             |
			 *       ____v0,1_____
			 *       |           |
			 *       |           |
			 * <- u0,0  H/P0,0   u1,0 ->
			 *       |           |
			 *   q0,0|___________|
			 *           v0,0
			 *             |
			 *             V
			 *
			 * V_t + q N x (P V) + grad( g P + 1/2 V*V) = 0
			 * P_t + div(P V) = 0
			 */
			/*
			 * U and V updates
			 */
			if(!param_nonlinear) //linear case
			{
				U = simVars.setup.h0*i_u;
				V = simVars.setup.h0*i_v;
			}
			else // nonlinear case
			{
			U = op.avg_b_x(i_h)*i_u;
			V = op.avg_b_y(i_h)*i_v;
			}

			if(!param_nonlinear) //linear case
				H = simVars.sim.g*i_h;// + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));
			else //non linear case
				H = simVars.sim.g*i_h + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));

			if(!param_nonlinear) //linear case
			{
				o_u_t = op.avg_f_y(simVars.sim.f0*op.avg_b_x(i_v)) - op.diff_b_x(H);
				o_v_t = -op.avg_f_x(simVars.sim.f0*op.avg_b_y(i_u)) - op.diff_b_y(H);
			}
			else //non linear case
			{
				q = (op.diff_b_x(i_v) - op.diff_b_y(i_u) + simVars.sim.f0) / op.avg_b_x(op.avg_b_y(i_h));

				o_u_t = op.avg_f_y(q*op.avg_b_x(V)) - op.diff_b_x(H);
				o_v_t = -op.avg_f_x(q*op.avg_b_y(U)) - op.diff_b_y(H);
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

//					double hx = simVars.disc.cell_size[0];
//					double hy = simVars.disc.cell_size[1];

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
				 * We use the hew v and u values to compute the update for p
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

			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
					o_ostream << "\tDIFF_P\tDIFF_U\tDIFF_V";

				if (param_compute_error){
					o_ostream << "\tANAL_DIFF_RMS_P\tANAL_DIFF_RMS_U\tANAL_DIFF_RMS_V";
					o_ostream << "\tANAL_DIFF_MAX_P\tANAL_DIFF_MAX_U\tANAL_DIFF_MAX_V";
				}
				o_ostream << std::endl;
			}

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

			o_ostream << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;

			if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
			{
				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// h
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
					}

				//benchmark_diff_h = (prog_h-tmp).reduce_norm1_quad() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
				benchmark_diff_h = (prog_h-tmp).reduce_maxAbs() ;
				o_ostream << "\t" << benchmark_diff_h;

				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// u space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
					}

				//benchmark_diff_u = (prog_u-tmp).reduce_norm1_quad() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
				benchmark_diff_u = (prog_u-tmp).reduce_maxAbs();
				o_ostream << "\t" << benchmark_diff_u;

				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
					}

				//benchmark_diff_v = (prog_v-tmp).reduce_norm1_quad() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
				benchmark_diff_v = (prog_v-tmp).reduce_maxAbs();
				o_ostream << "\t" << benchmark_diff_v;
			}

			if (param_compute_error)
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
		DataArray<2> t_h = t0_prog_h;
		DataArray<2> t_u = t0_prog_u;
		DataArray<2> t_v = t0_prog_v;

		rexiSWE.run_timestep_direct_solution(t_h, t_u, t_v, simVars.timecontrol.current_simulation_time, op, simVars);

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
	 * Arrays for online visualization and their textual description
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
				tmp = t_h;
				break;

			case -2:
				tmp = t_h-prog_h;	// difference
				break;
			}

			*o_dataArray = &tmp;
			*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
			return;
		}

#if 1
		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		tmp=**o_dataArray;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
#else
		DataArray<2> t_h = t0_prog_h;
		DataArray<2> t_u = t0_prog_u;
		DataArray<2> t_v = t0_prog_v;

		rexiSWE.run_timestep_direct_solution(t_h, t_u, t_v, simVars.timecontrol.current_simulation_time, op, simVars);

		switch(simVars.misc.vis_id)
		{
		case 0:
			tmp = prog_h;
			break;
		case 1:
			tmp = t_h;
			break;
		case 2:
			tmp = prog_u;
			break;
		case 3:
			tmp = t_u;
			break;
		case 4:
			tmp = prog_v;
			break;
		case 5:
			tmp = t_v;
			break;
		}


		*o_dataArray = &tmp;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
#endif
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
				description = "Direct solution for h";
				break;


			case -2:
				description = "Error in h";
				break;
			}
		}

		static char title_string[2048];
		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.6e, Vis: %s, Mass: %.6e, Energy: %.6e, Potential Enstrophy: %.6e, Max: %f, Min: %f ",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy,
				tmp.reduce_max(),
				tmp.reduce_min()	);

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
			"initial-freq-x-mul",
			"initial-freq-y-mul",
			"boundary-id",
			"rexi-zero-before-solving",
			"nonlinear",
			nullptr
	};


	simVars.bogus.var[0] = 0.2;
	simVars.bogus.var[1] = 256;	// M
	simVars.bogus.var[2] = 0;	// L = 0: default
	simVars.bogus.var[3] = 1;	// param_rexi_half
	simVars.bogus.var[4] = 0;
	simVars.bogus.var[5] = 0;
	simVars.bogus.var[6] = 0;
	simVars.bogus.var[7] = 1;	// Use spec diff per default for complex array
	simVars.bogus.var[8] = 0;	// Use spectral solver id
	simVars.bogus.var[9] = 1e-7;	// Error threshold

	simVars.bogus.var[10] = -1;	// freq. multiplier for initial conditions
	simVars.bogus.var[11] = -1;

	simVars.bogus.var[12] = 0;
	simVars.bogus.var[13] = 1;
	simVars.bogus.var[14] = 0;


	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--rexi-h [h-value]	h-sampling distance for REXI" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--rexi-m [m-value]	M-value for REXI (related to number of poles)" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--rexi-half [0/1]	Reduce rexi computations to its half" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "	--timestepping-mode [0/1/2]	Timestepping method to use" << std::endl;
		std::cout << "	                            0: Finite-difference / Spectral time stepping" << std::endl;
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
		std::cout << "	--rexi-helmholtz-solver-id=[int]	Use iterative solver for REXI" << std::endl;
		std::cout << "	--rexi-helmholtz-solver-eps=[err]	Error threshold for iterative solver" << std::endl;
		std::cout << std::endl;
		std::cout << "	--initial-freq-x-mul=[float]	Setup initial conditions by using this multiplier values" << std::endl;
		std::cout << "	--initial-freq-y-mul=[float]	" << std::endl;
		std::cout << std::endl;
		std::cout << "	--boundary-id=[0,1,...]	    Boundary id" << std::endl;
		std::cout << "                              0: no boundary" << std::endl;
		std::cout << "                              1: centered box" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-zero-before-solving=[0/1]	Zero the solution for the iterative solver" << std::endl;
		std::cout << std::endl;
		std::cout << "	--nonlinear=[0/1]	0-use linear equations, 1-use nonlinear eq" << std::endl;
		std::cout << std::endl;
		return -1;
	}

	param_rexi_h = simVars.bogus.var[0];
	param_rexi_m = simVars.bogus.var[1];
	param_rexi_l = simVars.bogus.var[2];
	param_rexi_half = simVars.bogus.var[3];
	param_timestepping_mode = simVars.bogus.var[4];
	param_compute_error = simVars.bogus.var[5];
	param_use_staggering = simVars.bogus.var[6];
	param_rexi_use_spectral_differences_for_complex_array = simVars.bogus.var[7];
	param_rexi_helmholtz_solver_id = simVars.bogus.var[8];
	param_rexi_helmholtz_solver_eps = simVars.bogus.var[9];

	param_initial_freq_x_mul = simVars.bogus.var[10];
	param_initial_freq_y_mul = simVars.bogus.var[11];

	param_boundary_id = simVars.bogus.var[12];

	param_rexi_zero_before_solving = simVars.bogus.var[13];

	param_nonlinear = simVars.bogus.var[14];

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
#if SWEET_GUI
		if (simVars.misc.gui_enabled)
		{
			VisSweet<SimulationSWE> visSweet(simulationSWE);
		}
		else
#endif

		{
			simulationSWE->reset();

			Stopwatch time;


			double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

			if (simVars.misc.verbosity > 1)
			{
				simulationSWE->update_diagnostics();
				diagnostics_energy_start = simVars.diag.total_energy;
				diagnostics_mass_start = simVars.diag.total_mass;
				diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
			}

#if SWEET_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			time.reset();


			while(true)
			{
				if (simVars.misc.verbosity > 1)
				{
					simulationSWE->timestep_output(buf);

					std::string output = buf.str();
					buf.str("");

					std::cout << output;

					if (simVars.misc.verbosity > 2)
						std::cerr << output;
				}

				if (simulationSWE->should_quit())
					break;

				simulationSWE->run_timestep();

				if (simulationSWE->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}
			}

			time.stop();

			double seconds = time();

			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;
			std::cout << "REXI alpha.size(): " << simulationSWE->rexiSWE.rexi.alpha.size() << std::endl;

			if (simVars.misc.verbosity > 1)
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
