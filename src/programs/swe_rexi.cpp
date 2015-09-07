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
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include "rexi/RexiSWE.hpp"

SimulationVariables simVars;

double param_rexi_h;
double param_rexi_m;
double param_rexi_l;
bool param_rexi_half;
int param_timestepping_mode;
bool param_compute_error;
bool param_use_staggering;
bool param_use_finite_differences_for_complex_array;


class SimulationSWE
{
public:
	DataArray<2> prog_h, prog_u, prog_v;
	DataArray<2> eta;
	DataArray<2> tmp;


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
		tmp(simVars.disc.res),

		t0_prog_h(simVars.disc.res),
		t0_prog_u(simVars.disc.res),
		t0_prog_v(simVars.disc.res),

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_diffs)
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

		prog_h.set_all(simVars.setup.h0);
		prog_u.set_all(0);
		prog_v.set_all(0);

		if (simVars.setup.scenario == -1)
		{
			std::cout << "Setting up steady state" << std::endl;

			if (simVars.sim.f == 0)
			{
				std::cerr << "Coriolis frequency f is set to 0" << std::endl;
				exit(1);
			}

			if (std::isinf(param_rexi_h))
				std::cerr << "ERROR: Steady state scenario here only makes sense for linear solver only" << std::endl;

			/**
			 * setup steady state
			 */
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
#if 0
					double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
#else
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
#endif

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
					double u = -simVars.sim.g/simVars.sim.f*bary;
					double v = simVars.sim.g/simVars.sim.f*barx;
#else
					double foo = std::exp(-50.0*(dx*dx + dy*dy));
					double h = simVars.setup.h0+foo;
					double u = -simVars.sim.g/simVars.sim.f*(-50.0*2.0*dy)*foo;
					double v = simVars.sim.g/simVars.sim.f*(-50.0*2.0*dx)*foo;
#endif

					prog_h.set(j, i, h);
					prog_u.set(j, i, u);
					prog_v.set(j, i, v);
				}
			}
		}
		else
		{
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

					prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
					prog_u.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
					prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
				}
			}
		}

		if (simVars.setup.input_data_filenames.size() > 0)
			prog_h.file_loadData_ascii(simVars.setup.input_data_filenames[0].c_str());

		if (simVars.setup.input_data_filenames.size() > 1)
			prog_u.file_loadData_ascii(simVars.setup.input_data_filenames[1].c_str());

		if (simVars.setup.input_data_filenames.size() > 2)
			prog_v.file_loadData_ascii(simVars.setup.input_data_filenames[2].c_str());

		t0_prog_h = prog_h;
		t0_prog_u = prog_u;
		t0_prog_v = prog_v;


		if (param_timestepping_mode == 1)
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
			rexiSWE.setup(-simVars.sim.CFL, param_rexi_h, param_rexi_m, param_rexi_l, simVars.sim.f, simVars.disc.res, simVars.sim.domain_size, param_rexi_half, param_use_finite_differences_for_complex_array);

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
		eta = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + simVars.sim.f) / prog_h;
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
			/*
			 * linearized non-conservative (advective) formulation:
			 *
			 * h_t = -h0*u_x - h0*v_y
			 * u_t = -g * h_x + f*v
			 * v_t = -g * h_y - f*u
			 */
			o_u_t = -simVars.sim.g*op.diff_c_x(i_h) + simVars.sim.f*i_v;
			o_v_t = -simVars.sim.g*op.diff_c_y(i_h) - simVars.sim.f*i_u;

			if (simVars.sim.viscosity != 0)
			{
				o_u_t -= op.diff2(i_u)*simVars.sim.viscosity;
				o_v_t -= op.diff2(i_v)*simVars.sim.viscosity;
			}

			if (simVars.sim.hyper_viscosity != 0)
			{
				o_u_t -= op.diff4(i_u)*simVars.sim.hyper_viscosity;
				o_v_t -= op.diff4(i_v)*simVars.sim.hyper_viscosity;
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
					double hx = simVars.disc.cell_size[0];
					double hy = simVars.disc.cell_size[1];

					double limit_speed = std::min(hx/i_u.reduce_maxAbs(), hy/i_v.reduce_maxAbs());

					// limit by re
					double limit_visc = std::numeric_limits<double>::infinity();
					if (simVars.sim.viscosity > 0)
						limit_visc = (hx*hx*hy*hy)/(4.0*simVars.sim.viscosity*simVars.sim.viscosity);
					if (simVars.sim.hyper_viscosity > 0)
						limit_visc = std::min((hx*hx*hx*hx*hy*hy*hy*hy)/(16.0*simVars.sim.hyper_viscosity*simVars.sim.hyper_viscosity), limit_visc);

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

					if (simVars.sim.viscosity != 0)
						o_h_t -= op.diff2(i_h)*simVars.sim.viscosity;

					if (simVars.sim.hyper_viscosity != 0)
						o_h_t -= op.diff2(i_h)*simVars.sim.hyper_viscosity;
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

			if (simVars.sim.potential_viscosity != 0)
				o_h_t -= op.diff2(i_h)*simVars.sim.potential_viscosity;

			if (simVars.sim.potential_hyper_viscosity != 0)
				o_h_t -= op.diff4(i_h)*simVars.sim.potential_hyper_viscosity;
		}
		else
		{
			DataArray<2> U(i_h.resolution), V(i_h.resolution), q(i_h.resolution), H(i_h.resolution);

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
//			U = op.avg_b_x(i_h)*i_u;
//			V = op.avg_b_y(i_h)*i_v;
			U = simVars.setup.h0*i_u;
			V = simVars.setup.h0*i_v;

			H = simVars.sim.g*i_h;// + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));

			o_u_t = op.avg_f_y(simVars.sim.f*op.avg_b_x(i_v)) - op.diff_b_x(H);
			o_v_t = -op.avg_f_x(simVars.sim.f*op.avg_b_y(i_u)) - op.diff_b_y(H);


			/*
			 * VISCOSITY
			 */
			if (simVars.sim.viscosity != 0)
			{
				o_u_t -= op.diff2(i_u)*simVars.sim.viscosity;
				o_v_t -= op.diff2(i_v)*simVars.sim.viscosity;
			}
			if (simVars.sim.hyper_viscosity != 0)
			{
				o_u_t -= op.diff4(i_u)*simVars.sim.hyper_viscosity;
				o_v_t -= op.diff4(i_v)*simVars.sim.hyper_viscosity;
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

					double hx = simVars.disc.cell_size[0];
					double hy = simVars.disc.cell_size[1];

					// limit by viscosity
					double limit_visc = std::numeric_limits<double>::infinity();
					if (simVars.sim.viscosity > 0)
						limit_visc = (hx*hx*hy*hy)/(4.0*simVars.sim.viscosity*simVars.sim.viscosity);
					if (simVars.sim.hyper_viscosity > 0)
						limit_visc = std::min((hx*hx*hx*hx*hy*hy*hy*hy)/(16.0*simVars.sim.hyper_viscosity*simVars.sim.hyper_viscosity), limit_visc);

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
				 * a kind of leapfrog:
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


			if (simVars.sim.potential_viscosity != 0)
				o_h_t -= op.diff2(i_h)*simVars.sim.potential_viscosity;

			if (simVars.sim.potential_hyper_viscosity != 0)
				o_h_t -= op.diff4(i_h)*simVars.sim.potential_hyper_viscosity;
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
			// REXI time stepping
			o_dt = -simVars.sim.CFL;
			rexiSWE.run_timestep(
					prog_h, prog_u, prog_v,
					op,
					simVars
			);
		}
		else if (param_timestepping_mode == 2)
		{
			// Analytical solution
			o_dt = -simVars.sim.CFL;
			rexiSWE.run_timestep_direct_solution(
					prog_h, prog_u, prog_v,
					-simVars.sim.CFL,
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
				o_ostream << "T\tMASS\tENERGY\tPOT_ENSTROPHY";

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
					o_ostream << "\tDIFF_P\tDIFF_U\tDIFF_V";

				if (param_compute_error)
					o_ostream << "\tANAL_DIFF_P\tANAL_DIFF_U\tANAL_DIFF_V";

				o_ostream << std::endl;
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

				benchmark_diff_h = (prog_h-tmp).reduce_norm1_quad() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
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

				benchmark_diff_u = (prog_u-tmp).reduce_norm1_quad() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
				o_ostream << "\t" << benchmark_diff_u;

				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));
					}

				benchmark_diff_v = (prog_v-tmp).reduce_norm1_quad() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
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

//		benchmark_analytical_error_maxabs_h = (t_h-prog_h).reduce_norm2_quad()/(double)(parameters.discretization.res[0]*parameters.discretization.res[1]);
//		benchmark_analytical_error_maxabs_u = (t_u-prog_u).reduce_norm2_quad()/(double)(parameters.discretization.res[0]*parameters.discretization.res[1]);
//		benchmark_analytical_error_maxabs_v = (t_v-prog_v).reduce_norm2_quad()/(double)(parameters.discretization.res[0]*parameters.discretization.res[1]);
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
#if 0
		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
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

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		static char title_string[2048];
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %.14s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				vis_arrays[id].description,
				simVars.diag.total_mass, simVars.diag.total_energy, simVars.diag.total_potential_enstrophy);

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

		case 'l':
			// dump data arrays
			prog_h.file_loadData_ascii("swe_rexi_dump_h.csv");
			prog_u.file_loadData_ascii("swe_rexi_dump_u.csv");
			prog_v.file_loadData_ascii("swe_rexi_dump_v.csv");
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
	const char *bogus_var_names[] = {
			"rexi-h",
			"rexi-m",
			"rexi-l",
			"rexi-half",
			"timestepping-mode",
			"compute-error",
			"staggering",
			"use-fd-for-complex-array",	/// use finite differences for complex array
			nullptr
	};


	simVars.bogus.var[0] = 0.2;
	simVars.bogus.var[1] = 256;	// M
	simVars.bogus.var[2] = 0;	// L = 0: default
	simVars.bogus.var[3] = 1;	// param_rexi_half
	simVars.bogus.var[4] = 0;
	simVars.bogus.var[5] = 0;
	simVars.bogus.var[6] = 0;
	simVars.bogus.var[7] = 0;	// don't use FD per default for complex array

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
		std::cout << "	--use-fd-for-complex-array=[0/1]	Use finite-differences for derivatives in spectral space" << std::endl;
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
	param_use_finite_differences_for_complex_array = simVars.bogus.var[7];


	SimulationSWE *simulationSWE = new SimulationSWE;

	std::ostringstream buf;
	buf << std::setprecision(14);

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
		time.reset();

		double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

		if (simVars.misc.verbosity > 1)
		{
			simulationSWE->update_diagnostics();
			diagnostics_energy_start = simVars.diag.total_energy;
			diagnostics_mass_start = simVars.diag.total_mass;
			diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
		}

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

	return 0;
}
