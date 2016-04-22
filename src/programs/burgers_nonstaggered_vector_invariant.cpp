
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

// AS tmp
#include <iostream>

SimulationVariables simVars;




double next_timestep_output = 0;


class SimulationSWECovariant
{
public:
	// prognostics
	DataArray<2> prog_u, prog_v;

	// temporary variables

	DataArray<2> tmp;

	Operators2D op;

	TimesteppingRK timestepping;

	int last_timestep_nr_update_diagnostics = -1;

	double benchmark_diff_u;
	double benchmark_diff_v;

	/**
	 * Two dimensional Burgers equation
	 *
	 * Equations:
	 *
	 *     \f$ u_t + uu_x + vu_y = \nu (u_xx + u_yy) \f$
	 *     \f$ v_t + uv_x + vv_y = \nu (v_xx + v_yy) \f$
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
		prog_u(simVars.disc.res),	// velocity (x-direction)
		prog_v(simVars.disc.res),	// velocity (y-direction)

		tmp(simVars.disc.res),

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
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

		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		simVars.reset();

		prog_u.set_all(0);
		prog_v.set_all(0);

		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{

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
		}

		if (simVars.setup.input_data_filenames.size() > 0)
			prog_u.file_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 1)
			prog_v.file_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);

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


		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res[0]*(double)simVars.disc.res[1]);

		// diagnostics_mass
		simVars.diag.total_mass = -1;

		// diagnostics_energy
		simVars.diag.total_energy =
			0.5*(
				(
					(prog_u*prog_u) +
					(prog_v*prog_v)
				)
			).reduce_sum_quad() * normalization;

		// potential enstropy
		simVars.diag.total_potential_enstrophy = -1;
	}


	void set_source( DataArray<2> &o_u_t )
	{
		double t = simVars.timecontrol.current_simulation_time;
		double tp = 2.0*M_PIl;

		/*
		 * f(t,x,y) = 2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)+2*PI*sin(2*PI*k*x)*cos(2*PI*k*x)*sin^2(2*PI*k*t)
		 *          - nu(-4*PI^2*k*sin(2*PI*k*x)*sin(2*PI*k*t))
		 * matching to:
		 * u(t,x,y) = 1/k * sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 57)
		{
			double k = 5.0;
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = tp * std::sin(tp*k*x)*std::cos(tp*k*t)
								  + tp*std::sin(tp*k*x)*std::sin(tp*k*t) * std::cos(tp*k*x)*std::sin(tp*k*t)
								  + simVars.sim.viscosity * (tp*tp*k* std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 2*2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)+4*2*PI*sin(2*PI*k*x)*cos(2*PI*k*x)*sin^2(2*PI*k*t)
		 *          - 2*nu(-4*PI^2*k*sin(2*PI*k*x)*sin(2*PI*k*t))
		 * matching to:
		 * u(t,x,y) = 2/k * sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 59)
		{
			double k = 1.0;
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = 2*tp * std::sin(tp*k*x)*std::cos(tp*k*t)
								  + 4*tp*std::sin(tp*k*x)*std::sin(tp*k*t) * std::cos(tp*k*x)*std::sin(tp*k*t)
								  + 2*simVars.sim.viscosity * (tp*tp*k* std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 2*PI*sin(2*PI*x)*cos(2*PI*t)+2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)
		 *			+ [sin(2*PI*x)*sin(2*PI*t)+1/k*sin(2*PI*k*x)*sin(2*PI*k*t)]
		 *			* [2*PI*cos(2*PI*x)*sin(2*PI*t)+2*PI*cos(2*PI*k*x)*sin(2*PI*k*t)]
		 *          - NU*[-4*PI*PI*sin(2*PI*x)*sin(2*PI*t)
		 *          - 4*PI*PI*k*sin(2*PI*k*x)*sin(2*PI*k*t)]
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)*sin(2*PI*t)+1/k*sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 58)
		{
			double k = 5.0;
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = tp*std::sin(tp*x)*std::cos(tp*t)+tp*std::sin(tp*k*x)*std::cos(tp*k*t)
					              + (std::sin(tp*x)*std::sin(tp*t)+1/k*std::sin(tp*k*x)*std::sin(tp*k*t))
					              * (tp*std::cos(tp*x)*std::sin(tp*t)+tp*std::cos(tp*k*x)*std::sin(tp*k*t))
								  - simVars.sim.viscosity*(-tp*tp*std::sin(tp*x)*std::sin(tp*t)
								  - tp*tp*k*std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 1
		 * matching to:
		 * u(t,x,y) = t
		 */
		if (simVars.setup.scenario == 51)
		{
			o_u_t.set_all(1.0);
		}

		/*
		 * f(t,x,y) = 2*t
		 * matching to:
		 * u(t,x,y) = t^2
		 */
		if (simVars.setup.scenario == 52)
		{
			o_u_t.set_all(2.0*t);
		}

		/*
		 * f(t,x,y) = 3*t^2
		 * matching to:
		 * u(t,x,y) = t^3
		 */
		if (simVars.setup.scenario == 53)
		{
			o_u_t.set_all(3.0*t*t);
		}

		/*
		 * f(t,x,y) = sin(2*PI*x) + t*sin(2*PI*x)*t*cos(2*PI*x)*2*PI - NU*(-4*PI*PI*t*sin(2*PI*x))
		 * matching to:
		 * u(t,x,y) = t*sin(2*PI*x)
		 */
		if (simVars.setup.scenario == 54)
		{
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = std::sin(tp*x)+t*std::sin(tp*x)*t*std::cos(tp*x)*tp
							      - simVars.sim.viscosity*(-tp*tp*t*std::sin(tp*x));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 2*PI*cos(2*PI*t)
		 * matching to:
		 * u(t,x,y) = sin(2*PI*t)
		 */
		if (simVars.setup.scenario == 55)
		{
			o_u_t.set_all(tp*std::cos(tp*t));
		}

		/*
		 * f(t,x,y) = 2*PI*cos(2*PI*k*t)
		 * matching to:
		 * u(t,x,y) = 1/k*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 56)
		{
			double k=5.0;
			o_u_t.set_all(tp*std::cos(tp*k*t));
		}
	}

	/**
	 * Compute derivative for time stepping and store it to
	 * u_t and v_t
	 */
	void p_run_euler_timestep_update(
			const DataArray<2> &i_u,	///< prognostic variables
			const DataArray<2> &i_v,	///< prognostic variables

			DataArray<2> &o_u_t,		///< time updates
			DataArray<2> &o_v_t,		///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		/* 2D Burgers equation
		 * u_t + u*u_x + v*u_y = nu*(u_xx + u_yy) [+ f(t,x,y)]
		 * v_t + u*v_x + v*v_y = nu*(v_xx + v_yy) [+ g(t,x,y)]
		 */

		/*
		 * u and v updates
		 */

		// Source-Term for manufactured solution
		double t = simVars.timecontrol.current_simulation_time;
		int inp;
		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{

				/*
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
					double tp = 2.0*M_PIl;
					double k = 1;

					/*
					 * f(t,x,y) = 2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)+2*PI*sin(2*PI*k*x)*cos(2*PI*k*x)*sin^2(2*PI*k*t)
					 *          - nu(-4*PI^2*k*sin(2*PI*k*x)*sin(2*PI*k*t))
					 * matching to:
					 * u(t,x,y) = 1/k * sin(2*PI*k*x)*sin(2*PI*k*t)
					 *
					double tmpvar = tp * std::sin(tp*k*x)*std::cos(tp*k*t)
					              + tp*std::sin(tp*k*x)*std::sin(tp*k*t) * std::cos(tp*k*x)*std::sin(tp*k*t)
					              + simVars.sim.viscosity * (tp*tp*k* std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
				*/

				/*
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
					double tp = 2.0*M_PIl;
					double k = 10.0;

					/*
					 * f(t,x,y) = 2*PI*[sin(2*PI*x)*cos(2*PI*t)+sin(2*PI*k*x)*cos(2*PI*k*t)]
					 *          + 2*PI*sin(2*PI*x)*sin(2*PI*t)*cos(2*PI*x)*sin(2*PI*t)
					 *          + 2*PI*sin(2*PI*x)*sin(2*PI*t)*cos(2*PI*k*x)*sin(2*PI*k*t)
					 *          + 2*PI/k*sin(2*PI*k*x)*sin(2*PI*k*t)*cos(2*PI*x)*sin(2*PI*t)
					 *          + 2*PI/k*sin(2*PI*k*x)*sin(2*PI*k*t)*cos(2*PI*k*x)*sin(2*PI*k*t)
					 *          + 4*PI*PI*NU*[sin(2*PI*x)*sin(2*PI*t)+k*sin(2*PI*k*x)*sin(2*PI*k*t)]
					 * matching to:
					 * u(t,x,y) = sin(2*PI*x)*sin(2*PI*t)+1/k*sin(2*PI*k*x)*sin(2*PI*k*t)
					 *
					double tmpvar = tp*(std::sin(tp*x)*std::cos(tp*t)+std::sin(tp*k*x)*std::cos(tp*k*x));
					tmpvar += tp*std::sin(tp*x)*std::sin(tp*t)*std::cos(tp*x)*std::sin(tp*t);
					tmpvar += tp*std::sin(tp*x)*std::sin(tp*t)*std::cos(tp*k*x)*std::sin(tp*k*t);
					tmpvar += tp/k*std::sin(tp*k*x)*std::sin(tp*k*t)*std::cos(tp*x)*std::sin(tp*t);
					tmpvar += tp/k*std::sin(tp*k*x)*std::sin(tp*k*t)*std::cos(tp*k*x)*std::sin(tp*k*t);
					tmpvar += tp*tp*simVars.sim.viscosity*(std::sin(tp*x)*std::sin(tp*t)+k*std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
				*/

			}
		}

		set_source(o_u_t);


		//necessary?
#if SWEET_USE_SPECTRAL_SPACE==1
		o_u_t.requestDataInSpectralSpace();
#endif

		/*
		 * reset to this line if no source term is used
		 * o_u_t = simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));
		 */
		o_u_t += simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));
		o_u_t -= i_u*op.diff_c_x(i_u) + i_v*op.diff_c_y(i_u);

		o_v_t = simVars.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));
		o_v_t -= i_u*op.diff_c_x(i_v) + i_v*op.diff_c_y(i_v);


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

//				if (simVars.misc.verbosity > 2)
//					std::cout << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

				o_dt = simVars.sim.CFL*std::min(limit_speed, limit_visc);
			}
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
				&SimulationSWECovariant::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_u, prog_v,
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


				prog_u.file_saveData_ascii((ss+"_u.csv").c_str());
				prog_v.file_saveData_ascii((ss+"_v.csv").c_str());

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

	VisStuff vis_arrays[2] =
	{
			{&prog_u,	"u"},
			{&prog_v,	"v"}
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
			if (simVars.misc.vis_id >= 2)
				simVars.misc.vis_id = 0;
			break;

		case 'V':
			simVars.misc.vis_id--;
			if (simVars.misc.vis_id < 0)
				simVars.misc.vis_id = 1;
			break;
		}
	}


	bool instability_detected()
	{
		return !(prog_u.reduce_all_finite() && prog_v.reduce_all_finite());
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
				std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
			}
		}
	}

	delete simulationSWE;

	return 0;
}
