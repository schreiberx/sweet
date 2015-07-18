
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationParameters.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Stopwatch.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

SimulationParameters parameters;


class SimulationSWECovariant
{
public:
	// prognostics
	DataArray<2> prog_P, prog_u, prog_v;

	// temporary variables
	DataArray<2> H, U, V;
	DataArray<2> q;

	// parameters
	DataArray<2> f;

	DataArray<2> tmp;

	Operators2D op;

	TimesteppingRK timestepping;

	int last_timestep_nr_update_diagnostics = -1;

	double benchmark_diff_h;
	double benchmark_diff_u;
	double benchmark_diff_v;

	/**
	 * See "The Dynamics of Finite-Difference Models of the Shallow-Water Equations", Robert Sadourny
	 *
	 * Prognostic:
	 *     V_t + q N x (P V) + grad(P + 0.5 V.V) = 0
	 *     P_t + div(P V) = 0
	 *
	 * Potential vorticity:
	 *     q = rot (V) / P
	 *
	 *   ____u0,1_____
	 *   |           |
	 *   |           |
	 * v0,0   P0,0   v1,0
	 *   |           |
	 *   |___u0,0____|
	 */
public:
	SimulationSWECovariant(
	)	:
		prog_P(parameters.res),	// density/pressure
		prog_u(parameters.res),	// velocity (x-direction)
		prog_v(parameters.res),	// velocity (y-direction)

		H(parameters.res),	//
		U(parameters.res),	// mass flux (x-direction)
		V(parameters.res),	// mass flux (y-direction)

		q(parameters.res),
		f(parameters.res),

		tmp(parameters.res),

		op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs)
	{
		reset();
	}

	~SimulationSWECovariant()
	{
	}

	void rot_coord(double angle, double &x, double &y)
	{
		angle *= 2.0*M_PI/360.0;
		double nx = std::cos(angle)*x - std::sin(angle)*y;
		double ny = std::sin(angle)*x + std::cos(angle)*y;
		x = nx;
		y = ny;
	}

	void rot_vector(double angle, double &x, double &y)
	{
		angle *= 2.0*M_PI/360.0;
		double nx = std::cos(angle)*x - std::sin(angle)*y;
		double ny = std::sin(angle)*x + std::cos(angle)*y;

		x = nx;
		y = ny;
	}


	void reset()
	{
		last_timestep_nr_update_diagnostics = -1;

		benchmark_diff_h = 0;
		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		parameters.reset();

		prog_P.setAll(parameters.setup_h0);
		prog_u.setAll(0);
		prog_v.setAll(0);

		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				{
					// h
					double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

					prog_P.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
				}

				{
					// u space
					double x = (((double)i)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

					prog_u.set(j,i, SWEValidationBenchmarks::return_u(parameters, x, y));
				}

				{
					// v space
					double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
					double y = (((double)j)/(double)parameters.res[1])*parameters.sim_domain_size[1];

					prog_v.set(j, i, SWEValidationBenchmarks::return_v(parameters, x, y));
				}
			}
		}


		timestep_output();
	}



	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == parameters.status_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = parameters.status_timestep_nr;

		if (parameters.use_f_array)
		{
			q = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + f) / prog_P;
		}
		else
		{
			q = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + parameters.sim_f) / prog_P;
		}

		double normalization = (parameters.sim_domain_size[0]*parameters.sim_domain_size[1]) /
								((double)parameters.res[0]*(double)parameters.res[1]);

		// diagnostics_mass
		parameters.diagnostics_mass = prog_P.reduce_sum_quad() * normalization;

		// diagnostics_energy
		parameters.diagnostics_energy = 0.5*(
				prog_P*prog_P +
				prog_P*(prog_u*prog_u) +
				prog_P*(prog_v*prog_v)
			).reduce_sum_quad() * normalization;

		// potential enstropy
		parameters.diagnostics_potential_entrophy = 0.5*(q*q*(prog_P)).reduce_sum_quad() * normalization;

	}



	void compute_upwinding_P_updates(
			const DataArray<2> &i_P,		///< prognostic variables (at T=tn)
			const DataArray<2> &i_u,		///< prognostic variables (at T=tn+dt)
			const DataArray<2> &i_v,		///< prognostic variables (at T=tn+dt)

			DataArray<2> &o_P_t	///< time updates (at T=tn+dt)
	)
	{
		//             |                       |                       |
		// --v---------|-----------v-----------|-----------v-----------|
		//   h-1       u0          h0          u1          h1          u2
		//

		// same a above, but formulated in a finite-difference style
		o_P_t =
			(
				(
					// u is positive
					op.shift_right(i_P)*i_u.return_value_if_positive()	// inflow
					-i_P*op.shift_left(i_u.return_value_if_positive())					// outflow

					// u is negative
					+(i_P*i_u.return_value_if_negative())	// outflow
					-op.shift_left(i_P*i_u.return_value_if_negative())		// inflow
				)*(1.0/parameters.sim_cell_size[0])	// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(i_P)*i_v.return_value_if_positive()		// inflow
					-i_P*op.shift_down(i_v.return_value_if_positive())					// outflow

					// v is negative
					+(i_P*i_v.return_value_if_negative())	// outflow
					-op.shift_down(i_P*i_v.return_value_if_negative())	// inflow
				)*(1.0/parameters.sim_cell_size[1])
			);
	}



	/**
	 * Compute derivative for time stepping and store it to
	 * P_t, u_t and v_t
	 */
	void p_run_euler_timestep_update(
			const DataArray<2> &i_h,	///< prognostic variables
			const DataArray<2> &i_u,	///< prognostic variables
			const DataArray<2> &i_v,	///< prognostic variables

			DataArray<2> &o_h_t,	///< time updates
			DataArray<2> &o_u_t,	///< time updates
			DataArray<2> &o_v_t,	///< time updates

			double &o_dt,			///< time step restriction
			double i_fixed_dt = 0		///< if this value is not equal to 0, use this time step size instead of computing one
	)
	{
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
		U = i_h*i_u;
		V = i_h*i_v;

		H = parameters.sim_g*i_h + 0.5*(i_u*i_u + i_v*i_v);

		if (parameters.setup_scenario != 5)
			q = (op.diff_c_x(i_v) - op.diff_c_y(i_u) + parameters.sim_f) / i_h;
		else
			q = (op.diff_c_x(i_v) - op.diff_c_y(i_u) + f) / i_h;

		o_u_t = q*V - op.diff_c_x(H);
		o_v_t = -q*U - op.diff_c_y(H);


		/*
		 * VISCOSITY
		 */
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


		/*
		 * P UPDATE
		 */
		if (!parameters.timestepping_leapfrog_like_update)
		{
			if (!parameters.timestepping_up_and_downwinding)
			{
				// standard update
				o_h_t = -op.diff_c_x(U) - op.diff_c_y(V);
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
			if (!parameters.timestepping_up_and_downwinding)
			{
				// recompute U and V
				U = i_h*(i_u+o_dt*o_u_t);
				V = i_h*(i_v+o_dt*o_v_t);

				// update based on new u and v values
				o_h_t = -op.diff_c_x(U) - op.diff_c_y(V);
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
		timestepping.run_rk_timestep(
				this,
				&SimulationSWECovariant::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_P, prog_u, prog_v,
				dt,
				parameters.timestepping_timestep_size,
				parameters.timestepping_runge_kutta_order
			);

		// provide information to parameters
		parameters.status_simulation_timestep_size = dt;
		parameters.status_simulation_time += dt;
		parameters.status_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif
	}

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
					o_ostream << "\tABS_P_DT\tABS_U_DT\tABS_V_DT";

				o_ostream << std::endl;
			}

			o_ostream << parameters.status_simulation_time << "\t" << parameters.diagnostics_mass << "\t" << parameters.diagnostics_energy << "\t" << parameters.diagnostics_potential_entrophy;

			// this should be zero for the steady state test
			if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
			{
				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// h
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
					}

				benchmark_diff_h = (prog_P-tmp).reduce_norm1() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_h;

				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// u space
						double x = (((double)i)/(double)parameters.res[0])*parameters.sim_domain_size[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_u(parameters, x, y));
					}

				benchmark_diff_u = (prog_u-tmp).reduce_norm1() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_v;

				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_size[0];
						double y = (((double)j)/(double)parameters.res[1])*parameters.sim_domain_size[1];

						tmp.set(j,i, SWEValidationBenchmarks::return_v(parameters, x, y));
					}

				benchmark_diff_v = (prog_v-tmp).reduce_norm1() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_v;
			}

			o_ostream << std::endl;
		}
	}



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
	void vis_post_frame_processing(int i_num_iterations)
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

	VisStuff vis_arrays[7] =
	{
			{&prog_P,	"P"},
			{&prog_u,	"u"},
			{&prog_v,	"v"},
			{&H,		"H"},
			{&q,		"q"},
			{&U,		"U"},
			{&V,		"V"}
	};

	void vis_get_vis_data_array(
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		int id = parameters.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		*o_aspect_ratio = parameters.sim_domain_size[1] / parameters.sim_domain_size[0];
	}



	const char* vis_get_status_string()
	{
		update_diagnostics();

		int id = parameters.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		static char title_string[1024];
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
			if (parameters.vis_id >= 7)
				parameters.vis_id = 6;
			break;

		case 'V':
			parameters.vis_id--;
			if (parameters.vis_id < 0)
				parameters.vis_id = 0;
			break;
		}
	}


	bool instability_detected()
	{
		return !(prog_P.reduce_all_finite() && prog_u.reduce_all_finite() && prog_v.reduce_all_finite());
	}
};




int main(int i_argc, char *i_argv[])
{
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);

	parameters.setup(i_argc, i_argv);

	SimulationSWECovariant *simulationSWE = new SimulationSWECovariant;

	std::ostringstream buf;
	buf << std::setprecision(14);



#if SWEET_GUI
	VisSweet<SimulationSWECovariant> visSweet(simulationSWE);
#else
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
		std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs(parameters.diagnostics_energy-diagnostics_energy_start) << std::endl;
		std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs(parameters.diagnostics_mass-diagnostics_mass_start) << std::endl;
		std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs(parameters.diagnostics_potential_entrophy-diagnostics_potential_entrophy_start) << std::endl;

		if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
		{
			std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark_diff_h << std::endl;
			std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
			std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
		}
	}
#endif

	delete simulationSWE;

	return 0;
}
