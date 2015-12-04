
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
#include <iomanip>
#include <stdio.h>



SimulationVariables simVars;

double next_timestep_output = 0;

class SimulationSWE
{
public:
	DataArray<2> prog_h, prog_u, prog_v;
	// beta plane
	DataArray<2> beta_plane;

	DataArray<2> eta;
	DataArray<2> tmp;


	Operators2D op;

	TimesteppingRK timestepping;

	int last_timestep_nr_update_diagnostics = -1;

	double benchmark_diff_h;
	double benchmark_diff_u;
	double benchmark_diff_v;

public:
	SimulationSWE()	:
		prog_h(simVars.disc.res),
		prog_u(simVars.disc.res),
		prog_v(simVars.disc.res),

		beta_plane(simVars.disc.res),

		eta(simVars.disc.res),
		tmp(simVars.disc.res),

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_diffs)
	{
		reset();
	}


	void reset()
	{
		next_timestep_output = 0;

		last_timestep_nr_update_diagnostics = -1;

		benchmark_diff_h = 0;
		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		simVars.timecontrol.current_timestep_nr = 0;
		simVars.timecontrol.current_simulation_time = 0;

		prog_h.set_all(simVars.setup.h0);
		prog_u.set_all(0);
		prog_v.set_all(0);

		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
				double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

				prog_h.set(j, i, SWEValidationBenchmarks::return_h(simVars, x, y));
				prog_u.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
				prog_v.set(j, i, SWEValidationBenchmarks::return_v(simVars, x, y));

				{
					// beta plane
					double y_beta = (((double)j+0.5)/(double)simVars.disc.res[1]);
					beta_plane.set(j, i, simVars.sim.f0+simVars.sim.beta*y_beta);
				}
			}
		}



		if (simVars.setup.input_data_filenames.size() > 0)
			prog_h.file_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

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
		if (simVars.sim.beta != 0.0)
			eta = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + simVars.sim.beta) / prog_h;
		else
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
					op.shift_right(i_h)*i_u.return_value_if_positive()	// inflow
					-i_h*op.shift_left(i_u.return_value_if_positive())					// outflow

					// u is negative
					+(i_h*i_u.return_value_if_negative())	// outflow
					-op.shift_left(i_h*i_u.return_value_if_negative())		// inflow
				)*(1.0/simVars.disc.cell_size[0])	// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(i_h)*i_v.return_value_if_positive()		// inflow
					-i_h*op.shift_down(i_v.return_value_if_positive())					// outflow

					// v is negative
					+(i_h*i_v.return_value_if_negative())	// outflow
					-op.shift_down(i_h*i_v.return_value_if_negative())	// inflow
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
		if (simVars.sim.top_bottom_zero_v_velocity)
		{
			/**
			 * See doc/beta_plane_for_swe/beta_plane_for_swe.lyx
			 * for further information.
			 *
			 * We use some tricks to simulate a slip condition (no outflow) at the top and bottom layer.
			 *
			 * Note, that this has to be done here, since the the RK time step updates may generate
			 * non-zero velocities at the top and bottom layer in the stages if only zeroing the components
			 * for the very 1st stage of RK.
			 */

			/*
			 * here, we first handle the VELOCITY UPDATES
			 *
			 * WARNING: This probably does not conserve the quantities
			 */
			((DataArray<2>&)i_v).set_row(0, 0);		// zero velocities
			((DataArray<2>&)i_v).set_row(-1, 0);

			((DataArray<2>&)i_u).copy_row(1, 0);	// slip boundaries
			((DataArray<2>&)i_u).copy_row(-2, -1);

			((DataArray<2>&)i_h).copy_row(1, 0);	// make smooth
			((DataArray<2>&)i_h).copy_row(-2, -1);
		}

		/*
		 * non-conservative (advective) formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
		 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
		 */
		o_u_t = -simVars.sim.g*op.diff_c_x(i_h) - i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
		o_v_t = -simVars.sim.g*op.diff_c_y(i_h) - i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

		if (simVars.sim.beta == 0.0)
		{
			o_u_t += simVars.sim.f0*i_v;
			o_v_t -= simVars.sim.f0*i_u;
		}
		else
		{
			o_u_t += beta_plane*i_v;
			o_v_t -= beta_plane*i_u;

			o_v_t.set_row(0, 0);
			o_v_t.set_row(-1, 0);
		}

#if 0
		if (simVars.sim.viscosity != 0)
		{
			o_u_t -= op.diffN_x(i_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			o_v_t -= op.diffN_y(i_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
		}
#endif


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
				o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
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
			std::cerr << "not supported!" << std::endl;

		}
	}



	void run_timestep()
	{
		double dt;

		// either set time step size to 0 for autodetection or to
		// a positive value to use a fixed time step size
		simVars.timecontrol.current_timestep_size = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 0);

		timestepping.run_rk_timestep(
				this,
				&SimulationSWE::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_h, prog_u, prog_v,
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



public:
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
				int secs = std::floor(simVars.timecontrol.current_simulation_time);
				int msecs = std::floor(1000000.*(simVars.timecontrol.current_simulation_time - floor(simVars.timecontrol.current_simulation_time)));
				char t_buf[256];
				sprintf(	t_buf,
							"%08d.%06d",
							secs, msecs
					);

				std::string ss = simVars.misc.output_file_name_prefix+"_t"+t_buf;

				prog_h.file_saveData_ascii((ss+"_h.csv").c_str());
				prog_u.file_saveData_ascii((ss+"_u.csv").c_str());
				prog_v.file_saveData_ascii((ss+"_v.csv").c_str());

				(op.diff_c_x(prog_v) - op.diff_c_y(prog_u)).file_saveData_ascii((ss+"_q.csv").c_str());
			}

			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
					o_ostream << "\tDIFF_P\tDIFF_U\tDIFF_V";

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

			o_ostream << std::endl;
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
		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
#if 0
		DataArray<2> &o = (DataArray<2> &)*vis_arrays[id].data;
		o.spec_setAll(0, 0);
		o.spec_set(0, 0, 1, 0);

		if (id == 0)
		{
			o.spec_set(0, 0, 1, 0);
		}
		else if (id == 1)
		{
			o.spec_set(0, 0, 1, 1);
		}
		else if (id == 2)
		{
			o.spec_set(0, 1, 1, 0);
		}
		else if (id == 3)
		{
			o.spec_set(0, 1, 1, 1);
//			o.spec_set(0, o.resolution_spec[0]-1, 1, 0);
		}
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
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
	{
		return -1;
	}


	SimulationSWE *simulationSWE = new SimulationSWE;

	std::ostringstream buf;

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

		while (true)
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
	}

	delete simulationSWE;

	return 0;
}
