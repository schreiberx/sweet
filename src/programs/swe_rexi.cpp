
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationParameters.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Stopwatch.hpp>
#include "../programs/RexiSWE.hpp"

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

SimulationParameters parameters;


double param_rexi_h;
double param_rexi_m;
bool param_rexi_half;
int param_timestepping_mode;
bool param_compute_error;


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

	double benchmark_analytical_diff_h;
	double benchmark_analytical_diff_u;
	double benchmark_analytical_diff_v;

	RexiSWE rexiSWE;


public:
	SimulationSWE()	:
		prog_h(parameters.res),
		prog_u(parameters.res),
		prog_v(parameters.res),

		eta(parameters.res),
		tmp(parameters.res),

		t0_prog_h(parameters.res),
		t0_prog_u(parameters.res),
		t0_prog_v(parameters.res),

		op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs)
	{
		reset();
	}


	void reset()
	{
		last_timestep_nr_update_diagnostics = -1;

		benchmark_diff_h = 0;
		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		benchmark_analytical_diff_h = 0;
		benchmark_analytical_diff_u = 0;
		benchmark_analytical_diff_v = 0;

		parameters.status_timestep_nr = 0;
		parameters.status_simulation_time = 0;

		prog_h.setAll(parameters.setup_h0);
		prog_u.setAll(0);
		prog_v.setAll(0);

		if (parameters.setup_scenario == -1)
		{
			std::cout << "Setting up steady state" << std::endl;

			if (parameters.sim_f == 0)
			{
				std::cerr << "Coriolis frequency f is set to 0" << std::endl;
				exit(1);
			}

			if (std::isinf(param_rexi_h))
			{
				std::cerr << "ERROR: Steady state scenario here only makes sense for linear solver only" << std::endl;
//				exit(1);
			}

			/**
			 * setup steady state
			 */
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
#if 0
					double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];
#else
					double x = (((double)i)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j)/(double)parameters.res[1])*parameters.sim_domain_size[1];
#endif

					// Gaussian
					double dx = x-parameters.setup_coord_x*parameters.sim_domain_size[0];
					double dy = y-parameters.setup_coord_y*parameters.sim_domain_size[1];

					double radius = parameters.setup_radius_scale*sqrt((double)parameters.sim_domain_size[0]*(double)parameters.sim_domain_size[0]+(double)parameters.sim_domain_size[1]*(double)parameters.sim_domain_size[1]);
					dx /= radius;
					dy /= radius;

#if 1
					double sx = x/parameters.sim_domain_size[0];
					double sy = y/parameters.sim_domain_size[1];

					double foo = std::sin(sx*2.0*M_PIl)*std::sin(sy*2.0*M_PIl);

					double barx = std::cos(sx*2.0*M_PIl)*std::sin(sy*2.0*M_PIl)*2.0*M_PIl/parameters.sim_domain_size[0];
					double bary = std::sin(sx*2.0*M_PIl)*std::cos(sy*2.0*M_PIl)*2.0*M_PIl/parameters.sim_domain_size[1];

					double h = parameters.setup_h0+foo;
					double u = -parameters.sim_g/parameters.sim_f*bary;
					double v = parameters.sim_g/parameters.sim_f*barx;
#else
					double foo = std::exp(-50.0*(dx*dx + dy*dy));
					double h = parameters.setup_h0+foo;
					double u = -parameters.sim_g/parameters.sim_f*(-50.0*2.0*dy)*foo;
					double v = parameters.sim_g/parameters.sim_f*(-50.0*2.0*dx)*foo;
#endif

					prog_h.set(j, i, h);
					prog_u.set(j, i, u);
					prog_v.set(j, i, v);
				}
			}
		}
		else
		{
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

					prog_h.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
					prog_u.set(j, i, SWEValidationBenchmarks::return_u(parameters, x, y));
					prog_v.set(j, i, SWEValidationBenchmarks::return_v(parameters, x, y));
				}
			}
		}

		t0_prog_h = prog_h;
		t0_prog_u = prog_u;
		t0_prog_v = prog_v;


		if (param_timestepping_mode == 1)
		{
			if (parameters.verbosity > 0)
			{
				std::cout << "REXI: Using REXI for time integration" << std::endl;
				std::cout << "REXI: Using " << param_rexi_m << " poles and sampling size " << param_rexi_h << std::endl;
			}

			if (parameters.sim_CFL >= 0)
			{
				std::cout << "Only constant time step size supported with REXI, use negative CFL to set constant time step size" << std::endl;
				exit(1);
			}

			// use REXI
			rexiSWE.setup(-parameters.sim_CFL, param_rexi_h, param_rexi_m, parameters.sim_f, parameters.res, parameters.sim_domain_size, param_rexi_half);

			if (parameters.verbosity > 2)
			{
				std::cout << "ALPHA:" << std::endl;
				for (std::size_t n = 0; n < rexiSWE.rexi.alpha.size(); n++)
					std::cout << rexiSWE.rexi.alpha[n] << std::endl;


				std::cout << "BETA:" << std::endl;
				for (std::size_t n = 0; n < rexiSWE.rexi.beta_re.size(); n++)
					std::cout << rexiSWE.rexi.beta_re[n] << std::endl;
			}
		}


		if (parameters.gui_enabled)
			timestep_output();
	}



	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == parameters.status_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = parameters.status_timestep_nr;

		double normalization = (parameters.sim_domain_size[0]*parameters.sim_domain_size[1]) /
								((double)parameters.res[0]*(double)parameters.res[1]);

		// mass
		parameters.diagnostics_mass = prog_h.reduce_sum_quad() * normalization;

		// energy
		parameters.diagnostics_energy = 0.5*(
				prog_h*prog_h +
				prog_h*prog_u*prog_u +
				prog_h*prog_v*prog_v
			).reduce_sum_quad() * normalization;

		// potential enstropy
		eta = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + parameters.sim_f) / prog_h;
		parameters.diagnostics_potential_entrophy = 0.5*(eta*eta*prog_h).reduce_sum_quad() * normalization;
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
					parameters.setup_h0*i_u.return_value_if_positive()	// inflow
					-parameters.setup_h0*op.shift_left(i_u.return_value_if_positive())					// outflow

					// u is negative
					+(parameters.setup_h0*i_u.return_value_if_negative())	// outflow
					-op.shift_left(parameters.setup_h0*i_u.return_value_if_negative())		// inflow
				)*(1.0/parameters.sim_cell_size[0])	// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					parameters.setup_h0*i_v.return_value_if_positive()		// inflow
					-parameters.setup_h0*op.shift_down(i_v.return_value_if_positive())					// outflow

					// v is negative
					+(parameters.setup_h0*i_v.return_value_if_negative())	// outflow
					-op.shift_down(parameters.setup_h0*i_v.return_value_if_negative())	// inflow
				)*(1.0/parameters.sim_cell_size[1])
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
		/*
		 * linearized non-conservative (advective) formulation:
		 *
		 * h_t = -h0*u_x - h0*v_y
		 * u_t = -g * h_x + f*v
		 * v_t = -g * h_y - f*u
		 */
		o_u_t = -parameters.sim_g*op.diff_c_x(i_h) + parameters.sim_f*i_v;
		o_v_t = -parameters.sim_g*op.diff_c_y(i_h) - parameters.sim_f*i_u;

		if (parameters.sim_viscocity != 0)
		{
			o_u_t += (op.diff2_c_x(i_u) + op.diff2_c_x(i_v))*parameters.sim_viscocity;
			o_v_t += (op.diff2_c_y(i_u) + op.diff2_c_y(i_v))*parameters.sim_viscocity;
		}

		if (parameters.sim_hyper_viscocity != 0)
		{
			o_u_t += (op.diff2_c_x(op.diff2_c_x(i_u)) + op.diff2_c_x(op.diff2_c_x(i_v)))*parameters.sim_hyper_viscocity;
			o_v_t += (op.diff2_c_y(op.diff2_c_y(i_u)) + op.diff2_c_y(op.diff2_c_y(i_v)))*parameters.sim_hyper_viscocity;
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
				double limit_speed = std::min(parameters.sim_cell_size[0]/i_u.reduce_maxAbs(), parameters.sim_cell_size[1]/i_v.reduce_maxAbs());

				double hx = parameters.sim_cell_size[0];
				double hy = parameters.sim_cell_size[1];

				// limit by re
				double limit_visc = std::numeric_limits<double>::infinity();
				if (parameters.sim_viscocity > 0)
					limit_visc = (hx*hx*hy*hy)/(4.0*parameters.sim_viscocity*parameters.sim_viscocity);
				if (parameters.sim_hyper_viscocity > 0)
					limit_visc = std::min((hx*hx*hx*hx*hy*hy*hy*hy)/(16.0*parameters.sim_hyper_viscocity*parameters.sim_hyper_viscocity), limit_visc);

				// limit by gravitational acceleration
				double limit_gh = std::min(parameters.sim_cell_size[0], parameters.sim_cell_size[1])/std::sqrt(parameters.sim_g*i_h.reduce_maxAbs());

				if (parameters.verbosity > 2)
					std::cerr << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

				o_dt = parameters.sim_CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);
			}
		}

		if (!parameters.timestepping_leapfrog_like_update)
		{
			if (!parameters.timestepping_up_and_downwinding)
			{
				// standard update
				o_h_t = -op.diff_c_x(i_u)*parameters.setup_h0 - op.diff_c_y(i_v*parameters.setup_h0);
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

				if (parameters.sim_viscocity != 0)
					o_h_t += (op.diff2_c_x(i_h) + op.diff2_c_y(i_h))*parameters.sim_viscocity;

				if (parameters.sim_hyper_viscocity != 0)
					o_h_t += (op.diff2_c_x(op.diff2_c_x(i_h)) + op.diff2_c_y(op.diff2_c_y(i_h)))*parameters.sim_hyper_viscocity;
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
			if (!parameters.timestepping_up_and_downwinding)
			{
				// recompute U and V

				// update based on new u and v values
				o_h_t = -op.diff_c_x(
							parameters.setup_h0*(i_u+o_dt*o_u_t)
						)
						- op.diff_c_y(
								parameters.setup_h0*(i_v+o_dt*o_v_t)
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

		if (parameters.sim_potential_viscocity != 0)
			o_h_t += (op.diff2_c_x(i_h) + op.diff2_c_y(i_h))*parameters.sim_potential_viscocity;

		if (parameters.sim_potential_hyper_viscocity != 0)
			o_h_t += (op.diff2_c_x(op.diff2_c_x(i_h)) + op.diff2_c_y(op.diff2_c_y(i_h)))*parameters.sim_potential_hyper_viscocity;

	}



	void run_timestep()
	{
		double dt;
		if (param_timestepping_mode == 0)
		{
			if (parameters.sim_CFL < 0)
				parameters.timestepping_timestep_size = -parameters.sim_CFL;

			// standard time stepping
			timestepping.run_rk_timestep(
					this,
					&SimulationSWE::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					prog_h, prog_u, prog_v,
					dt,
					parameters.timestepping_timestep_size,
					parameters.timestepping_runge_kutta_order,
					parameters.status_simulation_time
				);
		}
		else if (param_timestepping_mode == 1)
		{
			// REXI time stepping
			dt = -parameters.sim_CFL;
			rexiSWE.run_timestep(
					prog_h, prog_u, prog_v,
					op,
					parameters
			);
		}
		else if (param_timestepping_mode == 2)
		{
			// Analytical solution
			dt = -parameters.sim_CFL;
			rexiSWE.run_timestep_direct_solution(
					prog_h, prog_u, prog_v,
					-parameters.sim_CFL,
					op,
					parameters
			);
		}
		else
		{
			std::cerr << "Invalid time stepping method" << std::endl;
			exit(1);
		}

		// provide information to parameters
		parameters.status_simulation_timestep_size = dt;
		parameters.status_simulation_time += dt;
		parameters.status_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif
	}



public:
	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (parameters.verbosity > 0)
		{
			update_diagnostics();

			if (parameters.status_timestep_nr == 0)
			{
				o_ostream << "T\tMASS\tENERGY\tPOT_ENSTROPHY";

				if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
					o_ostream << "\tDIFF_P\tDIFF_U\tDIFF_V";

				if (param_compute_error)
					o_ostream << "\tANAL_DIFF_P\tANAL_DIFF_U\tANAL_DIFF_V";

				o_ostream << std::endl;
			}

			o_ostream << parameters.status_simulation_time << "\t" << parameters.diagnostics_mass << "\t" << parameters.diagnostics_energy << "\t" << parameters.diagnostics_potential_entrophy;

			if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
			{
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// h
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
					}

				benchmark_diff_h = (prog_h-tmp).reduce_norm1_quad() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_h;

				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// u space
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_u(parameters, x, y));
					}

				benchmark_diff_u = (prog_u-tmp).reduce_norm1_quad() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_u;

				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_v(parameters, x, y));
					}

				benchmark_diff_v = (prog_v-tmp).reduce_norm1_quad() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_v;
			}

			if (param_compute_error)
			{
				DataArray<2> t_h = t0_prog_h;
				DataArray<2> t_u = t0_prog_u;
				DataArray<2> t_v = t0_prog_v;

				rexiSWE.run_timestep_direct_solution(t_h, t_u, t_v, parameters.status_simulation_time, op, parameters);

				benchmark_analytical_diff_h = (t_h-prog_h).reduce_norm2_quad() / (double)(parameters.res[0]*parameters.res[1]);
				benchmark_analytical_diff_u = (t_u-prog_u).reduce_norm2_quad() / (double)(parameters.res[0]*parameters.res[1]);
				benchmark_analytical_diff_v = (t_v-prog_v).reduce_norm2_quad() / (double)(parameters.res[0]*parameters.res[1]);

				o_ostream << "\t" << benchmark_analytical_diff_h;
				o_ostream << "\t" << benchmark_analytical_diff_u;
				o_ostream << "\t" << benchmark_analytical_diff_v;
			}

			o_ostream << std::endl;
		}
	}



public:
	bool should_quit()
	{
		if (parameters.max_timesteps_nr != -1 && parameters.max_timesteps_nr <= parameters.status_timestep_nr)
			return true;

		if (parameters.max_simulation_time != -1 && parameters.max_simulation_time <= parameters.status_simulation_time)
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
		if (parameters.run_simulation)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}


	struct VisStuff
	{
		const DataArray<2>* data;
		const char *description;
	};

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
		int id = parameters.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		*o_aspect_ratio = parameters.sim_domain_size[1] / parameters.sim_domain_size[0];
#else
		DataArray<2> t_h = t0_prog_h;
		DataArray<2> t_u = t0_prog_u;
		DataArray<2> t_v = t0_prog_v;

		rexiSWE.run_timestep_direct_solution(t_h, t_u, t_v, parameters.status_simulation_time, op, parameters);

		switch(parameters.vis_id)
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
		*o_aspect_ratio = parameters.sim_domain_size[1] / parameters.sim_domain_size[0];
#endif
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostic values if required
		update_diagnostics();

		int id = parameters.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		static char title_string[2048];
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %.14s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				parameters.status_simulation_time,
				parameters.status_simulation_time/(60.0*60.0*24.0),
				parameters.status_timestep_nr,
				parameters.status_simulation_timestep_size,
				vis_arrays[id].description,
				parameters.diagnostics_mass, parameters.diagnostics_energy, parameters.diagnostics_potential_entrophy);

		return title_string;
	}



	void vis_pause()
	{
		parameters.run_simulation = !parameters.run_simulation;
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			parameters.vis_id++;
			break;

		case 'V':
			parameters.vis_id--;
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
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);

	const char *bogus_var_names[] = {
			"rexi-h",
			"rexi-m",
			"rexi-half",
			"timestepping-mode",
			"compute-error",
			nullptr
	};


	parameters.bogus_var0 = 0.1;
	parameters.bogus_var1 = 200;
	parameters.bogus_var2 = 1;	// param_rexi_half
	parameters.bogus_var3 = 0;
	parameters.bogus_var4 = 0;

	if (!parameters.setup(i_argc, i_argv, bogus_var_names))
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
		return -1;
	}

	param_rexi_h = parameters.bogus_var0;
	param_rexi_m = parameters.bogus_var1;
	param_rexi_half = parameters.bogus_var2;
	param_timestepping_mode = parameters.bogus_var3;
	param_compute_error = parameters.bogus_var4;

	SimulationSWE *simulationSWE = new SimulationSWE;

	std::ostringstream buf;
	buf << std::setprecision(14);

#if SWEET_GUI
	if (parameters.gui_enabled)
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

		if (parameters.verbosity > 1)
		{
			simulationSWE->update_diagnostics();
			diagnostics_energy_start = parameters.diagnostics_energy;
			diagnostics_mass_start = parameters.diagnostics_mass;
			diagnostics_potential_entrophy_start = parameters.diagnostics_potential_entrophy;
		}

		while(true)
		{
			if (parameters.verbosity > 1)
			{
				simulationSWE->timestep_output(buf);

				std::string output = buf.str();
				buf.str("");

				std::cout << output;

				if (parameters.verbosity > 2)
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

		std::cout << "Simulation time: " << seconds << " seconds" << std::endl;
		std::cout << "Time per time step: " << seconds/(double)parameters.status_timestep_nr << " sec/ts" << std::endl;

		if (parameters.verbosity > 1)
		{
			std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((parameters.diagnostics_energy-diagnostics_energy_start)/diagnostics_energy_start) << std::endl;
			std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((parameters.diagnostics_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
			std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((parameters.diagnostics_potential_entrophy-diagnostics_potential_entrophy_start)/diagnostics_potential_entrophy_start) << std::endl;

			if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
			{
				std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark_diff_h << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
			}

			if (param_compute_error)
			{
				std::cout << "DIAGNOSTICS ANALYTICAL DIFF H:\t" << simulationSWE->benchmark_analytical_diff_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL DIFF U:\t" << simulationSWE->benchmark_analytical_diff_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL DIFF V:\t" << simulationSWE->benchmark_analytical_diff_v << std::endl;
			}
		}
	}

	delete simulationSWE;

	return 0;
}
