
#include "../include/sweet/plane/PlaneData.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDiagnostics.hpp>
#include <sweet/Stopwatch.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;


SimulationVariables simVars;

double next_timestep_output = 0;


class SimulationSWECovariant
{
public:
	// prognostics
	PlaneData prog_h, prog_u, prog_v;

	// temporary variables
	PlaneData H, U, V;
	PlaneData q;			/// vorticity

	// beta plane
	PlaneData beta_plane;

	PlaneData tmp;

	PlaneOperators op;

	PlaneDataTimesteppingRK timestepping;

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
	 *   |            |
	 *   v0,0 P0,0 v1,0
	 *   |            |
	 *   |    u0,0    |
	 *   |____________|
	 */
public:
	SimulationSWECovariant(
	)	:
		prog_h(planeDataConfig),	// density/pressure
		prog_u(planeDataConfig),	// velocity (x-direction)
		prog_v(planeDataConfig),	// velocity (y-direction)

		H(planeDataConfig),	//
		U(planeDataConfig),	// mass flux (x-direction)
		V(planeDataConfig),	// mass flux (y-direction)

		q(planeDataConfig),
		beta_plane(planeDataConfig),

		tmp(planeDataConfig),

		op(
				planeDataConfig,
				simVars.sim.domain_size,
				simVars.disc.use_spectral_basis_diffs
		)
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

		prog_h.physical_set_all(simVars.sim.h0);
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);




		prog_h.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
				io_data = SWEPlaneBenchmarks::return_h(simVars, x, y);
			}
		);

		prog_u.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double x = (((double)i)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
				io_data = SWEPlaneBenchmarks::return_u(simVars, x, y);
			}
		);

		prog_v.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
				double y = (((double)j)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
				io_data = SWEPlaneBenchmarks::return_v(simVars, x, y);
			}
		);
/*
		beta_plane.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				double y_beta = (((double)j+0.5)/(double)simVars.disc.res_physical[1]);
				io_data = simVars.sim.f0+simVars.sim.beta*y_beta;
			}
		);
*/
		if (simVars.setup.input_data_filenames.size() > 0)
			prog_h.file_physical_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 1)
			prog_u.file_physical_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 2)
			prog_v.file_physical_loadData(simVars.setup.input_data_filenames[2].c_str(), simVars.setup.input_data_binary);

		if (simVars.misc.gui_enabled)
			timestep_output();
	}



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



	void compute_upwinding_P_updates(
			const PlaneData &i_P,		///< prognostic variables (at T=tn)
			const PlaneData &i_u,		///< prognostic variables (at T=tn+dt)
			const PlaneData &i_v,		///< prognostic variables (at T=tn+dt)

			PlaneData &o_P_t	///< time updates (at T=tn+dt)
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
					op.shift_right(i_P)*i_u.physical_query_return_value_if_positive()	// inflow
					-i_P*op.shift_left(i_u.physical_query_return_value_if_positive())					// outflow

					// u is negative
					+(i_P*i_u.physical_query_return_value_if_negative())	// outflow
					-op.shift_left(i_P*i_u.physical_query_return_value_if_negative())		// inflow
				)*(1.0/simVars.disc.cell_size[0])	// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(i_P)*i_v.physical_query_return_value_if_positive()		// inflow
					-i_P*op.shift_down(i_v.physical_query_return_value_if_positive())					// outflow

					// v is negative
					+(i_P*i_v.physical_query_return_value_if_negative())	// outflow
					-op.shift_down(i_P*i_v.physical_query_return_value_if_negative())	// inflow
				)*(1.0/simVars.disc.cell_size[1])
			);
	}



	/**
	 * Compute derivative for time stepping and store it to
	 * P_t, u_t and v_t
	 */
	void p_run_euler_timestep_update(
			const PlaneData &i_h,	///< prognostic variables
			const PlaneData &i_u,	///< prognostic variables
			const PlaneData &i_v,	///< prognostic variables

			PlaneData &o_h_t,		///< time updates
			PlaneData &o_u_t,		///< time updates
			PlaneData &o_v_t,		///< time updates

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

		H = simVars.sim.gravitation*i_h + 0.5*(i_u*i_u + i_v*i_v);

		if (simVars.setup.benchmark_scenario_id != 5)
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

//				double hx = simVars.disc.cell_size[0];
//				double hy = simVars.disc.cell_size[1];

				// limit by re
				double limit_visc = std::numeric_limits<double>::infinity();

				// limit by gravitational acceleration
				double limit_gh = std::min(simVars.disc.cell_size[0], simVars.disc.cell_size[1])/std::sqrt(simVars.sim.gravitation*i_h.reduce_maxAbs());

//				if (simVars.misc.verbosity > 2)
//					std::cout << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

				o_dt = simVars.sim.CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);
			}
		}

		/*
		 * P UPDATE
		 */
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

		timestepping.run_timestep(
				this,
				&SimulationSWECovariant::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_h, prog_u, prog_v,
				dt,
				simVars.timecontrol.current_timestep_size,
				simVars.disc.timestepping_order,
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
				write_file(prog_h, "prog_P");
				write_file(prog_u, "prog_u");
				write_file(prog_v, "prog_v");

				write_file(op.diff_c_x(prog_v) - op.diff_c_y(prog_u), "prog_q");
			}

			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				if (simVars.setup.benchmark_scenario_id == 2 || simVars.setup.benchmark_scenario_id == 3 || simVars.setup.benchmark_scenario_id == 4)
					o_ostream << "\tABS_P_DT\tABS_U_DT\tABS_V_DT";

				o_ostream << std::endl;
			}

			o_ostream << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;

			// this should be zero for the steady state test
			if (simVars.setup.benchmark_scenario_id == 2 || simVars.setup.benchmark_scenario_id == 3 || simVars.setup.benchmark_scenario_id == 4)
			{

				tmp.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
					{
						double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

						io_data = SWEPlaneBenchmarks::return_h(simVars, x, y);
					}
				);

				benchmark_diff_h = (prog_h-tmp).reduce_norm1_quad() / (double)(simVars.disc.res_physical[0]*simVars.disc.res_physical[1]);
				o_ostream << "\t" << benchmark_diff_h;


				// set data to something to overcome assertion error

				tmp.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
					{
						// u space
						double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

						io_data = SWEPlaneBenchmarks::return_u(simVars, x, y);
					}
				);

				benchmark_diff_u = (prog_u-tmp).reduce_norm1_quad() / (double)(simVars.disc.res_physical[0]*simVars.disc.res_physical[1]);
				o_ostream << "\t" << benchmark_diff_u;


				tmp.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

						io_data = SWEPlaneBenchmarks::return_v(simVars, x, y);
					}
				);

				benchmark_diff_v = (prog_v-tmp).reduce_norm1() / (double)(simVars.disc.res_physical[0]*simVars.disc.res_physical[1]);
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
		const PlaneData* data;
		const char *description;
	};

	VisStuff vis_arrays[7] =
	{
			{&prog_h,	"P"},
			{&prog_u,	"u"},
			{&prog_v,	"v"},
			{&H,		"H"},
			{&q,		"q"},
			{&U,		"U"},
			{&V,		"V"}
	};

	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data
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
		sprintf(title_string, "Time: %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %.14s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
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
		return !(prog_h.reduce_boolean_all_finite() && prog_u.reduce_boolean_all_finite() && prog_v.reduce_boolean_all_finite());
	}
};




int main(int i_argc, char *i_argv[])
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

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

			if (simVars.setup.benchmark_scenario_id == 2 || simVars.setup.benchmark_scenario_id == 3 || simVars.setup.benchmark_scenario_id == 4)
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
