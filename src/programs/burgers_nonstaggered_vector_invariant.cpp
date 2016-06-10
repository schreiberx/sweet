
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

#include <parareal/Parareal.hpp>

SimulationVariables simVars;




double next_timestep_output = 0;
bool param_imex = 1;


class SimulationInstance
#if SWEET_PARAREAL
		:
		public Parareal_SimulationInstance
#endif
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
	SimulationInstance(
	)	:
		prog_u(simVars.disc.res),	// velocity (x-direction)
		prog_v(simVars.disc.res),	// velocity (y-direction)

		tmp(simVars.disc.res),

		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
#if SWEET_PARAREAL != 0
		,
		_parareal_data_start_u(simVars.disc.res), _parareal_data_start_v(simVars.disc.res),
		_parareal_data_fine_u(simVars.disc.res), _parareal_data_fine_v(simVars.disc.res),
		_parareal_data_coarse_u(simVars.disc.res), _parareal_data_coarse_v(simVars.disc.res),
		_parareal_data_output_u(simVars.disc.res), _parareal_data_output_v(simVars.disc.res),
		_parareal_data_error_u(simVars.disc.res), _parareal_data_error_v(simVars.disc.res)
#endif
	{
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
			double k = simVars.sim.f0;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
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
			double k = simVars.sim.f0;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
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
			double k = simVars.sim.f0;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
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
		 * f(t,x,y) = 1000*sin(2*PI*x) + 1000^2*t*sin(2*PI*x)*t*cos(2*PI*x)*2*PI - 1000*NU*(-4*PI*PI*t*sin(2*PI*x))
		 * matching to:
		 * u(t,x,y) = 1000*t*sin(2*PI*x)
		 */
		if (simVars.setup.scenario == 54)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = 1000*std::sin(tp*x)+1000*1000*t*std::sin(tp*x)*t*std::cos(tp*x)*tp
							      - 1000*simVars.sim.viscosity*(-tp*tp*t*std::sin(tp*x));

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
			double k=simVars.sim.f0;
			o_u_t.set_all(tp*std::cos(tp*k*t));
		}

		/*
		 * f(t,x,y) = sin(2*PI*x)*cos(2*PI*x)*2*PI - NU*(-4*PI*PI*sin(2*PI*x))
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)
		 */
		if (simVars.setup.scenario == 60)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = std::sin(tp*x)*std::cos(tp*x)*tp
							      - simVars.sim.viscosity*(-tp*tp*std::sin(tp*x));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/* Test for 2u-grad(u)=F
		 * f(t,x,y) = (1+4*PI*PI)sin(2*PI*x)
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)
		 */
		if (simVars.setup.scenario == 61)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = (1+tp*tp)*std::sin(tp*x);

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 0.25*[SUM (cos(**)*(-PI*k)+sin(**)*(2*PI*k)^2)*# + 0.25*SUM sin(**)*# * SUM cos(**)*2*PI*k*#]
		 * mit: ** = 2*PI*k*x-PI*k*t+PI*k
		 * und: #  = EPS/sinh(0.5*PI*k*EPS)
		 * matching to:
		 * u(t,x,y) = 0.5*SUM_(k=1)^(k_max) sin(2*PI*k*x-PI*k*t+PI*k)*EPS/sinh(0.5*PI*k*EPS)
		 */
		if (simVars.setup.scenario == 62)
		{
			int kmax = 100;
			double eps = 0.1;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double tmpvar = 0;
					for (int k = 1; k < kmax; k++)
					{
						double argument = tp*k*x - M_PIl*k*t + M_PIl*k;
						tmpvar += cos(argument)*(-M_PIl*k) + sin(argument)*(tp*k)*(tp*k)*eps/sinh(M_PIl*k*eps/2);
						for (int l = 1; l < kmax; l++)
						{
							double argument2 = tp*l*x - M_PIl*l*t + M_PIl*l;
							tmpvar += sin(argument)*eps/sinh(M_PIl*k*eps/2)*cos(argument2)*eps/sinh(M_PIl*l*eps/2)*tp*l;
						}
					}

					tmpvar *= 0.25;
					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		if (simVars.setup.scenario == 63)
			o_u_t.set_all(0.0);

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
		/* 2D Burgers equation [with source term]
		 * u_t + u*u_x + v*u_y = nu*(u_xx + u_yy) [+ f(t,x,y)]
		 * v_t + u*v_x + v*v_y = nu*(v_xx + v_yy) [+ g(t,x,y)]
		 */

		DataArray<2> f(i_u.resolution);
		set_source(f);

		/*
		 * u and v updates
		 */

		o_u_t = -(i_u*op.diff_c_x(i_u) + i_v*op.diff_c_y(i_u));
		o_u_t += simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));
		// Delete this line if no source is used.
		o_u_t += f;

		o_v_t = -(i_u*op.diff_c_x(i_v) + i_v*op.diff_c_y(i_v));
		o_v_t += simVars.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));

		/*
		 * TIME STEP SIZE
		 */
		if (i_fixed_dt > 0)
		{
			o_dt = i_fixed_dt;
		}
		else
		{
			std::cout << "Only fixed time step size supported" << std::endl;
			assert(false);
			exit(1);
		}

	}



	void run_timestep()
	{
		double dt;

		// either set time step size to 0 for autodetection or to
		// a positive value to use a fixed time step size
		simVars.timecontrol.current_timestep_size = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 0);

		if (!param_imex)
		{
			// run standard RK
			timestepping.run_rk_timestep(
					this,
					&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					prog_u, prog_v,
					dt,
					simVars.timecontrol.current_timestep_size,
					simVars.disc.timestepping_runge_kutta_order,
					simVars.timecontrol.current_simulation_time
				);
		}
		else
		{
			run_timestep_imex(
					prog_u, prog_v,
					simVars.timecontrol.current_timestep_size,
					op,
					simVars
			);
		}

		dt = simVars.timecontrol.current_timestep_size;

		// provide information to parameters
		simVars.timecontrol.current_timestep_size = dt;
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

//#if SWEET_GUI
		timestep_output();
//#endif
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
				// not necessary for 1D Problems
				//prog_v.file_saveData_ascii((ss+"_v.csv").c_str());

				// set data to something to overcome assertion error
				/*for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// u space
						double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
					}

				(prog_u-tmp).file_saveData_ascii((ss+"_u_diff.csv").c_str());
				*/

				/* Not necessary for 1D problems
				for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						tmp.set(j,i, SWEValidationBenchmarks::return_v(simVars, x, y));
					}

				(prog_v-tmp).file_saveData_ascii((ss+"_v_diff.csv").c_str());
				*/

			}

			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
					o_ostream << "\tABS_P_DT\tABS_U_DT\tABS_V_DT";

				o_ostream << std::endl;

			}

			o_ostream << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;

		}
		o_ostream << std::endl;
	}


public:
	void compute_errors()
	{
		/*TODO: rewrite to compute errors between analytical and actual solution
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
		*/
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


#if SWEET_PARAREAL

	/******************************************************
	 ******************************************************
	 *       ************** PARAREAL **************
	 ******************************************************
	 ******************************************************/

	DataArray<2> _parareal_data_start_u, _parareal_data_start_v;
	Parareal_Data_DataArrays<2> parareal_data_start;

	DataArray<2> _parareal_data_fine_u, _parareal_data_fine_v;
	Parareal_Data_DataArrays<2> parareal_data_fine;

	DataArray<2> _parareal_data_coarse_u, _parareal_data_coarse_v;
	Parareal_Data_DataArrays<2> parareal_data_coarse;

	DataArray<2> _parareal_data_output_u, _parareal_data_output_v;
	Parareal_Data_DataArrays<2> parareal_data_output;

	DataArray<2> _parareal_data_error_u, _parareal_data_error_v;
	Parareal_Data_DataArrays<2> parareal_data_error;

	double timeframe_start = -1;
	double timeframe_end = -1;

	bool output_data_valid = false;

	void parareal_setup()
	{
		{
			DataArray<2>* data_array[2] = {&_parareal_data_start_u, &_parareal_data_start_v};
			parareal_data_start.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_fine_u, &_parareal_data_fine_v};
			parareal_data_fine.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_coarse_u, &_parareal_data_coarse_v};
			parareal_data_coarse.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_output_u, &_parareal_data_output_v};
			parareal_data_output.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_error_u, &_parareal_data_error_v};
			parareal_data_error.setup(data_array);
		}

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


		*parareal_data_start.data_arrays[0] = prog_u;
		*parareal_data_start.data_arrays[1] = prog_v;

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

		prog_u = *parareal_data_start.data_arrays[0];
		prog_v = *parareal_data_start.data_arrays[1];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			this->run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// copy to buffers
		*parareal_data_fine.data_arrays[0] = prog_u;
		*parareal_data_fine.data_arrays[1] = prog_v;
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

#endif

	/**
	 * IMEX time stepping for the coarse timestepping method
	 *
	 * This routine solves the system
	 * (I-\nu\Delta) u(t+\tau) = u(t) - \tau (u(t)*\nabla(u(t))
	 * for u(t+\tau).
	 */
	bool run_timestep_imex(
			DataArray<2> &io_u,
			DataArray<2> &io_v,

			double i_timestep_size,	///< timestep size

			Operators2D &op,
			const SimulationVariables &i_simVars
			)
	{
		DataArray<2> u=io_u;
		DataArray<2> v=io_v;

		// Initialize and set source for manufactured solution
		DataArray<2> f(io_u.resolution);
		set_source(f);

		// Modify timestep to final time if necessary
		double t = 0.0;
		if (simVars.timecontrol.current_simulation_time+i_timestep_size<simVars.timecontrol.max_simulation_time)
			t = i_timestep_size;
		else
			t = simVars.timecontrol.max_simulation_time-simVars.timecontrol.current_simulation_time;

		// Setting explicit right hand side and operator of the left hand side
		DataArray<2> rhs_u = u - t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + t*f;
		DataArray<2> rhs_v = v - t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));
		DataArray<2>   lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).addScalar_Cart(1.0);

#if 1   // solving the system directly by inverting the left hand side operator
		io_u = rhs_u.spec_div_element_wise(lhs);
		io_v = rhs_v.spec_div_element_wise(lhs);

#else	// making the second step of the IMEX-RK1 scheme
		DataArray<2> u1 = rhs_u.spec_div_element_wise(lhs);
		DataArray<2> v1 = rhs_v.spec_div_element_wise(lhs);

		io_u = u + t*simVars.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
				- t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) +f*t;
		io_v = v + t*simVars.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
				- t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));
#endif

		return true;
	}

#if SWEET_PARAREAL

	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	void run_timestep_coarse()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_coarse()" << std::endl;

		prog_u = *parareal_data_start.data_arrays[0];
		prog_v = *parareal_data_start.data_arrays[1];

		// run implicit time step
//		assert(i_max_simulation_time < 0);
//		assert(simVars.sim.CFL < 0);
		if (!param_imex)
		{
			// run standard RK
			double dt = 0.0;
			timestepping.run_rk_timestep(
					this,
					&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					prog_u, prog_v,
					dt,
					timeframe_end - timeframe_start,
					simVars.disc.timestepping_runge_kutta_order,
					simVars.timecontrol.current_simulation_time
				);
		}
		else
		{
			run_timestep_imex(
						prog_u, prog_v,
						timeframe_end - timeframe_start,
						op,
						simVars
				);
		}

		std::cout << "Run timestep coarse: " << prog_u.reduce_rms() << std::endl;

		/*
		// Test for lowpass filter
		if (prog_u.array_data_spectral_space_valid)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
			for (std::size_t i = prog_u.array_data_spectral_length/2; i < prog_u.array_data_spectral_length; i+=2)
			{
				double factor = std::sqrt(1/(1+i-prog_u.array_data_spectral_length/2));
				prog_u.array_data_spectral_space[i] *= factor;
				prog_u.array_data_spectral_space[i+1] *= factor;
			}
		}

		if (prog_v.array_data_spectral_space_valid)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
			for (std::size_t i = prog_u.array_data_spectral_length/2; i < prog_v.array_data_spectral_length; i+=2)
			{
				double factor = std::sqrt(1/(1+i-prog_u.array_data_spectral_length/2));
				prog_v.array_data_spectral_space[i] *= factor;
				prog_v.array_data_spectral_space[i+1] *= factor;
			}
		}
*/

		// copy to buffers
		*parareal_data_coarse.data_arrays[0] = prog_u;
		*parareal_data_coarse.data_arrays[1] = prog_v;
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

		for (int k = 0; k < 2; k++)
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
			for (int k = 0; k < 2; k++)
				*parareal_data_output.data_arrays[k] = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			output_data_valid = true;
			return convergence;
		}



		for (int k = 0; k < 2; k++)
		{
			tmp = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			convergence = std::max(
					convergence,
					(*parareal_data_output.data_arrays[k]-tmp).reduce_maxAbs()
				);

			*parareal_data_output.data_arrays[k] = tmp;
		}

		simVars.timecontrol.current_simulation_time = timeframe_end;
		prog_u = *parareal_data_output.data_arrays[0];
		prog_v = *parareal_data_output.data_arrays[1];
		compute_errors(); //TODO still to be implemented

		//std::cout << "maxabs error compared to analytical solution: " << benchmark_analytical_error_maxabs_h << std::endl;

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
		Parareal_Data_DataArrays<2>& data = (Parareal_Data_DataArrays<2>&)i_data;

		std::ostringstream ss;
		ss << simVars.misc.output_file_name_prefix << "output_iter" << iteration_id << "_slice" << time_slice_id << ".csv";

		std::string filename = ss.str();

//		std::cout << "filename: " << filename << std::endl;
		data.data_arrays[0]->file_saveData_ascii(filename.c_str());

		// set data to something to overcome assertion error
#if SWEET_THREADING
#pragma omp parallel for OPENMP_SIMD
#endif
		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
				// u space
				double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
				double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

				tmp.set(j, i, SWEValidationBenchmarks::return_u(simVars, x, y));
			}

		ss.str("");
		ss.clear();
		ss << simVars.misc.output_file_name_prefix << "output_iter" << iteration_id << "_slice" << time_slice_id << "_diff.csv";
		filename = ss.str();
		(*data.data_arrays[0]-tmp).file_saveData_ascii(filename.c_str());

		//data.data_arrays[0]->file_saveData_vtk(filename.c_str(), filename.c_str());
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
	NUMABlockAlloc::setup();

	const char *bogus_var_names[] = {
			"use-imex",
			nullptr
	};

	simVars.bogus.var[0] = param_imex;


	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--use-imex [0/1]	Use IMEX time stepping" << std::endl;
		std::cout << "" << std::endl;

		return -1;
	}
	param_imex = simVars.bogus.var[0];


	std::ostringstream buf;
	buf << std::setprecision(14);

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

#if SWEET_GUI
	if (simVars.misc.gui_enabled)
	{
		SimulationInstance *simulationBurgers = new SimulationInstance;
		VisSweet<SimulationInstance> visSweet(simulationBurgers);
		delete simulationBurgers;
	}
	else
#endif
	{
		SimulationInstance *simulationBurgers = new SimulationInstance;
		simulationBurgers->reset();

		Stopwatch time;
		time.reset();


		double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

		if (simVars.misc.verbosity > 1)
		{
			simulationBurgers->update_diagnostics();
			diagnostics_energy_start = simVars.diag.total_energy;
			diagnostics_mass_start = simVars.diag.total_mass;
			diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
		}

		while(true)
		{
			if (simVars.misc.verbosity > 1)
			{
				simulationBurgers->timestep_output(buf);

				std::string output = buf.str();
				buf.str("");

				std::cout << output << std::flush;
			}

			if (simulationBurgers->should_quit())
				break;

			simulationBurgers->run_timestep();

			if (simulationBurgers->instability_detected())
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
				std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationBurgers->benchmark_diff_u << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationBurgers->benchmark_diff_v << std::endl;
			}
		}

		delete simulationBurgers;
	}


	return 0;
}
