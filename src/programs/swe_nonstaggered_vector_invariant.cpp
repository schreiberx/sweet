
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Stopwatch.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

SimulationVariables simVars;

double next_timestep_output = 0;


class SimulationSWECovariant
{
public:
	// prognostics
	DataArray<2> prog_P, prog_u, prog_v;

	// temporary variables
	DataArray<2> H, U, V;
	DataArray<2> q;			/// vorticity

	// beta plane
	DataArray<2> beta_plane;

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
	 *
	 *     \f$ V_t + q N x (P V) + grad(P + 0.5 V.V) = 0 \f$
	 *
	 *     \f$ P_t + div(P V) = 0 \f$
	 *
	 * Potential vorticity:
	 *     \f$  q = rot (V) / P \f$
	 *
	 *   ______________
	 *   |            |
	 *   |    u0,1    |
	 *   v0,0 P0,0 v1,0
	 *   |    u0,0    |
	 *   |____________|
	 */
public:
	SimulationSWECovariant(
	)	:
		prog_P(simVars.disc.res),	// density/pressure
		prog_u(simVars.disc.res),	// velocity (x-direction)
		prog_v(simVars.disc.res),	// velocity (y-direction)

		H(simVars.disc.res),	//
		U(simVars.disc.res),	// mass flux (x-direction)
		V(simVars.disc.res),	// mass flux (y-direction)

		q(simVars.disc.res),
		beta_plane(simVars.disc.res),

		tmp(simVars.disc.res),

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_diffs)
	{
		reset();
	}


	~SimulationSWECovariant()
	{
	}


	void reset()
	{
		next_timestep_output = 0;

		last_timestep_nr_update_diagnostics = -1;

		benchmark_diff_h = 0;
		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		simVars.reset();

		prog_P.set_all(simVars.setup.h0);
		prog_u.set_all(0);
		prog_v.set_all(0);

		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
				{
					// h
					double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

					prog_P.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
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

				{
					// beta plane
					double y_beta = (((double)j+0.5)/(double)simVars.disc.res[1]);
					beta_plane.set(j, i, simVars.sim.f0+simVars.sim.beta*y_beta);
				}
			}
		}

		if (simVars.setup.input_data_filenames.size() > 0)
			prog_P.file_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 1)
			prog_u.file_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 2)
			prog_v.file_loadData(simVars.setup.input_data_filenames[2].c_str(), simVars.setup.input_data_binary);

		if (simVars.misc.gui_enabled)
			timestep_output();
	}



	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;

		if (simVars.sim.beta != 0)
		{
			q = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + beta_plane) / prog_P;
		}
		else
		{
			q = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + simVars.sim.f0) / prog_P;
		}

		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res[0]*(double)simVars.disc.res[1]);

		// diagnostics_mass
		simVars.diag.total_mass = prog_P.reduce_sum_quad() * normalization;

		// diagnostics_energy
		simVars.diag.total_energy =
			0.5*(
				prog_P*
				(
					prog_P +
					(prog_u*prog_u) +
					(prog_v*prog_v)
				)
			).reduce_sum_quad() * normalization;

		// potential enstropy
		simVars.diag.total_potential_enstrophy = 0.5*(q*q*(prog_P)).reduce_sum_quad() * normalization;
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
				)*(1.0/simVars.disc.cell_size[0])	// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(i_P)*i_v.return_value_if_positive()		// inflow
					-i_P*op.shift_down(i_v.return_value_if_positive())					// outflow

					// v is negative
					+(i_P*i_v.return_value_if_negative())	// outflow
					-op.shift_down(i_P*i_v.return_value_if_negative())	// inflow
				)*(1.0/simVars.disc.cell_size[1])
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

			DataArray<2> &o_h_t,		///< time updates
			DataArray<2> &o_u_t,		///< time updates
			DataArray<2> &o_v_t,		///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
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

		H = simVars.sim.g*i_h + 0.5*(i_u*i_u + i_v*i_v);

		if (simVars.setup.scenario != 5)
			q = (op.diff_c_x(i_v) - op.diff_c_y(i_u) + simVars.sim.f0) / i_h;
		else
			q = (op.diff_c_x(i_v) - op.diff_c_y(i_u) + beta_plane) / i_h;

		o_u_t = q*V - op.diff_c_x(H);
		o_v_t = -q*U - op.diff_c_y(H);


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

				double hx = simVars.disc.cell_size[0];
				double hy = simVars.disc.cell_size[1];

				// limit by re
				double limit_visc = std::numeric_limits<double>::infinity();

				// limit by gravitational acceleration
				double limit_gh = std::min(simVars.disc.cell_size[0], simVars.disc.cell_size[1])/std::sqrt(simVars.sim.g*i_h.reduce_maxAbs());

//				if (simVars.misc.verbosity > 2)
//					std::cout << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

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
			if (!simVars.disc.timestepping_up_and_downwinding)
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
#if 0
		if (simVars.sim.potential_viscosity != 0)
			o_h_t -= op.diff2(i_h)*simVars.sim.potential_viscosity;

		if (simVars.sim.potential_hyper_viscosity != 0)
			o_h_t -= op.diff4(i_h)*simVars.sim.potential_hyper_viscosity;
#endif
	}



	void run_timestep()
	{
		double dt;

		// either set time step size to 0 for autodetection or to
		// a positive value to use a fixed time step size
		simVars.timecontrol.current_timestep_size = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 0);

		timestepping.run_rk_timestep(
				this,
				&SimulationSWECovariant::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_P, prog_u, prog_v,
				dt,
				simVars.timecontrol.current_timestep_size,
				simVars.disc.timestepping_runge_kutta_order,
				simVars.timecontrol.current_simulation_time
			);

		// provide information to parameters
		simVars.timecontrol.current_timestep_size = dt;
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif
	}

	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (simVars.timecontrol.current_simulation_time < next_timestep_output)
			return;

		if (simVars.misc.be_verbose_after_this_simulation_time_period != 0)
		{
			// advance to next time step output
			while (next_timestep_output <= simVars.timecontrol.current_simulation_time)
				next_timestep_output += simVars.misc.be_verbose_after_this_simulation_time_period;
		}

		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();

			if (simVars.misc.output_file_name_prefix.size() > 0)
			{
				double secs = simVars.timecontrol.current_simulation_time;
				double msecs = 1000000.*(simVars.timecontrol.current_simulation_time - floor(simVars.timecontrol.current_simulation_time));
				char t_buf[256];
				sprintf(	t_buf,
							"%08d.%06d",
							(int)secs, (int)msecs
					);

				std::string ss = simVars.misc.output_file_name_prefix+"_t"+t_buf;


				prog_P.file_saveData_ascii((ss+"_h.csv").c_str());
				prog_u.file_saveData_ascii((ss+"_u.csv").c_str());
				prog_v.file_saveData_ascii((ss+"_v.csv").c_str());

				(op.diff_c_x(prog_v) - op.diff_c_y(prog_u)).file_saveData_ascii((ss+"_q.csv").c_str());
			}

			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
					o_ostream << "\tABS_P_DT\tABS_U_DT\tABS_V_DT";

				o_ostream << std::endl;
			}

			o_ostream << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;

			// this should be zero for the steady state test
			if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
			{
				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// h
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
					}

				benchmark_diff_h = (prog_P-tmp).reduce_norm1() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
				o_ostream << "\t" << benchmark_diff_h;

				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// u space
						double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
					}

				benchmark_diff_u = (prog_u-tmp).reduce_norm1() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
				o_ostream << "\t" << benchmark_diff_v;

				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j,i, SWEValidationBenchmarks::return_v(simVars, x, y));
					}

				benchmark_diff_v = (prog_v-tmp).reduce_norm1() / (double)(simVars.disc.res[0]*simVars.disc.res[1]);
				o_ostream << "\t" << benchmark_diff_v;
			}

			o_ostream << std::endl;
		}
	}



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
	void vis_post_frame_processing(int i_num_iterations)
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
		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
	}



	const char* vis_get_status_string()
	{
		update_diagnostics();

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		static char title_string[1024];
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
			if (simVars.misc.vis_id >= 7)
				simVars.misc.vis_id = 6;
			break;

		case 'V':
			simVars.misc.vis_id--;
			if (simVars.misc.vis_id < 0)
				simVars.misc.vis_id = 0;
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
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	SimulationSWECovariant *simulationSWE = new SimulationSWECovariant;

	std::ostringstream buf;
	buf << std::setprecision(14);


#if SWEET_GUI
	if (simVars.misc.gui_enabled)
	{
		VisSweet<SimulationSWECovariant> visSweet(simulationSWE);
	}
	else
#endif
	{
//		simulationSWE->reset();

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

				std::cout << output << std::flush;
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
		std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
		std::cout << "Timesteps: " << simVars.timecontrol.current_timestep_nr << std::endl;


		if (simVars.misc.verbosity > 1)
		{
			std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs(simVars.diag.total_energy-diagnostics_energy_start) << std::endl;
			std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs(simVars.diag.total_mass-diagnostics_mass_start) << std::endl;
			std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs(simVars.diag.total_potential_enstrophy-diagnostics_potential_entrophy_start) << std::endl;

			if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
			{
				std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark_diff_h << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
			}
		}
	}

	delete simulationSWE;

	return 0;
}
