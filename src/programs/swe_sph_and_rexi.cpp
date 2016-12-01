/*
 * AppTestSWE.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
	#include <sweet/plane/PlaneDataConfig.hpp>
	#include <sweet/plane/PlaneData.hpp>
	#include <sweet/Convert_SphereData_To_PlaneData.hpp>
#endif

#include <benchmarks_sphere/BenchmarkCombined.hpp>

#include <sweet/sphere/SphereData.hpp>

// explicit time stepping
#include <sweet/sphere/SphereDataTimesteppingExplicitRK.hpp>

// implicit time stepping for SWE
#include <sweet/sphere/app_swe/SWEImplicit_SPHRobert.hpp>

#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalComplex.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>


#include "swe_sphere_rexi/SWE_Sphere_REXI.hpp"



SimulationVariables simVars;

//Diagnostic measures at initial stage
//double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

// Plane data config
SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

// Plane data config
SphereDataConfig sphereDataConfigInstanceExt;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigInstanceExt;



#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif

/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */
bool param_rexi_use_coriolis_formulation = true;
bool param_compute_error = false;
double param_geostr_balance_freq_multiplier = 1.0;



class SimulationInstance
{
public:
	SphereOperators op;
	SphereOperatorsComplex opComplex;

	// Runge-Kutta stuff
	SphereDataTimesteppingExplicitRK timestepping_explicit;

	// Implicit timestepping solver
	SWEImplicit_SPHRobert timestepping_implicit_swe;

	SWE_Sphere_REXI swe_sphere_rexi;

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	SphereData prog_h;
	SphereData prog_u;
	SphereData prog_v;

	REXI rexi;

#if SWEET_GUI
	PlaneData viz_plane_data;
#endif

	int render_primitive_id;



public:
	SimulationInstance()	:
		prog_h(sphereDataConfig),
		prog_u(sphereDataConfig),
		prog_v(sphereDataConfig)

#if SWEET_GUI
		,viz_plane_data(planeDataConfig)
#endif
	{
		reset();
	}





	SphereData f(SphereData i_sphData)
	{
		return op.mu(i_sphData*2.0*simVars.sim.coriolis_omega);
	}


	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;

		// TODO: Calculate accurate normalization
		FatalError("TODO: calculate accurate normalization");
		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res_physical[0]*(double)simVars.disc.res_physical[1]);


		// mass
		simVars.diag.total_mass = prog_h.physical_reduce_sum_metric() * normalization;

		// energy
		simVars.diag.total_energy = 0.5*((
				prog_h*prog_h +
				prog_h*prog_u*prog_u +
				prog_h*prog_v*prog_v
			).physical_reduce_sum_metric()) * normalization;

		SphereData eta(sphereDataConfig);

		SphereData f(sphereDataConfig);
		double two_omega = 2.0*simVars.sim.coriolis_omega;
		f.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = mu*two_omega;
				}
			);

		// potential vorticity and pot. enstropy
		if (simVars.misc.sphere_use_robert_functions)
		{
			eta = (op.robert_vort(prog_u, prog_v) + f) / prog_h;
		}
		else
		{
			eta = (op.vort(prog_u, prog_v) + f) / prog_h;
		}

		simVars.diag.total_potential_enstrophy = 0.5*(eta*eta*prog_h).physical_reduce_sum_metric() * normalization;

	}


	void reset()
	{
		render_primitive_id = 1;

		// one month runtime
		if (simVars.timecontrol.max_simulation_time == -1)
		{
			simVars.timecontrol.max_simulation_time = 31*60*60*24;

			// 144 h
			//simVars.timecontrol.max_simulation_time = 144*60*60;

			// 200 h
			simVars.timecontrol.max_simulation_time = 200*60*60;
	//		simVars.timecontrol.max_simulation_time = 1;
		}

		// Diagnostics measures
		last_timestep_nr_update_diagnostics = -1;

		simVars.misc.output_next_sim_seconds = 0;


		if (simVars.sim.CFL < 0)
			simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;


		if (simVars.timecontrol.current_timestep_size <= 0)
		{
			// TRY to guess optimal time step size

			// time step size
			if (sphereDataConfig->physical_num_lat < 256)
				simVars.timecontrol.current_timestep_size = 0.002*simVars.sim.earth_radius/(double)sphereDataConfig->physical_num_lat;
			else
				simVars.timecontrol.current_timestep_size = 0.001*simVars.sim.earth_radius/(double)sphereDataConfig->physical_num_lat;
		}


		if (simVars.misc.output_each_sim_seconds <= 0)
		{
			simVars.misc.output_each_sim_seconds = 60*30;	// output every 1/2 hour
		}


		if (simVars.rexi.use_rexi == 1)
		{
			// Override for REXI
			if (simVars.timecontrol.current_timestep_size <= 0)
			{
				std::cout << "Timestep size not positive" << std::endl;
				assert(false);
				exit(1);
			}

			simVars.misc.use_nonlinear_equations = 0;
		}


		BenchmarkCombined::setupInitialConditions(prog_h, prog_u, prog_v, simVars, op);

		if (simVars.sim.coriolis_omega != 0)
			param_rexi_use_coriolis_formulation = true;

		std::cout << "Using time step size dt = " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << "Running simulation until t_end = " << simVars.timecontrol.max_simulation_time << std::endl;
		std::cout << "Parameters:" << std::endl;
		std::cout << " + Gravity: " << simVars.sim.gravitation << std::endl;
		std::cout << " + Earth_radius: " << simVars.sim.earth_radius << std::endl;
		std::cout << " + Average height: " << simVars.sim.h0 << std::endl;
		std::cout << " + Coriolis_omega: " << simVars.sim.coriolis_omega << std::endl;
		std::cout << " + Viscosity D: " << simVars.sim.viscosity << std::endl;
		std::cout << " + use_nonlinear: " << simVars.misc.use_nonlinear_equations << std::endl;
		std::cout << " + Use REXI Coriolis formulation: " << (param_rexi_use_coriolis_formulation ? "true" : "false") << std::endl;
		std::cout << std::endl;
		std::cout << " + Benchmark scenario id: " << simVars.setup.benchmark_scenario_id << std::endl;
		std::cout << " + Use robert functions: " << simVars.misc.sphere_use_robert_functions << std::endl;
		std::cout << " + Use REXI: " << simVars.rexi.use_rexi << " (0: explicit TS, 1: REXI, 2: implicit)" << std::endl;
		std::cout << " + REXI h: " << simVars.rexi.rexi_h << std::endl;
		std::cout << " + REXI M: " << simVars.rexi.rexi_M << std::endl;
		std::cout << " + REXI use half poles: " << simVars.rexi.rexi_use_half_poles << std::endl;
		std::cout << " + REXI normalization: " << simVars.rexi.rexi_normalization << std::endl;
		std::cout << " + REXI additional modes: " << simVars.rexi.rexi_use_extended_modes << std::endl;
		std::cout << std::endl;
		std::cout << " + RK order: " << simVars.disc.timestepping_runge_kutta_order << std::endl;
		std::cout << " + timestep size: " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << " + output timestep size: " << simVars.misc.output_each_sim_seconds << std::endl;

		std::cout << std::endl;

		if (simVars.rexi.use_rexi == 1)
		{
			swe_sphere_rexi.setup(
					simVars.rexi.rexi_h,
					simVars.rexi.rexi_M,
					simVars.rexi.rexi_L,

					sphereDataConfig,
					&simVars.sim,
					simVars.timecontrol.current_timestep_size,

					simVars.rexi.rexi_use_half_poles,
					simVars.misc.sphere_use_robert_functions,
					simVars.rexi.rexi_use_extended_modes,
					simVars.rexi.rexi_normalization,
					param_rexi_use_coriolis_formulation
				);
		}
		else if (simVars.rexi.use_rexi == -1)
		{
			if (simVars.sim.CFL >= 0)
				FatalError("CFL >= 0: Set negative CFL for constant time step size");

			timestepping_implicit_swe.setup(
					sphereDataConfig,
					sphereDataConfigExt,

					simVars.sim.earth_radius,
					simVars.sim.coriolis_omega,
					simVars.sim.gravitation*simVars.sim.h0,
					-simVars.sim.CFL,

					param_rexi_use_coriolis_formulation
				);
		}

	}



	/**
	 *
	 */
	std::string write_file(
			const SphereData &i_sphereData,
			const char* i_name,
			bool i_phi_shifted
		)
	{
		char buffer[1024];

		const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
		if (i_phi_shifted)
			i_sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			i_sphereData.physical_file_write(buffer);

		return buffer;
	}


	void write_file_output()
	{
		if (simVars.misc.output_file_name_prefix.length() == 0)
			return;

		std::string output_filename;

		std::cout << "Simulation time: " << simVars.timecontrol.current_simulation_time << std::endl;


		output_filename = write_file(prog_h, "h", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << " (min: " << prog_h.physical_reduce_min() << ", max: " << prog_h.physical_reduce_max() << ")" << std::endl;

		output_filename = write_file(prog_h, "u", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << std::endl;

		output_filename = write_file(prog_h, "v", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << std::endl;

		output_filename = write_file(op.vort(prog_u, prog_v)/simVars.sim.earth_radius, "eta", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << std::endl;


		/*
		 * write reference solution
		 */
		if (	param_compute_error &&
				simVars.misc.use_nonlinear_equations == 0 &&
				simVars.setup.benchmark_scenario_id == 4
		)
		{
			SphereData test_h(sphereDataConfig);
			SphereData test_u(sphereDataConfig);
			SphereData test_v(sphereDataConfig);

			BenchmarkCombined::setupInitialConditions(test_h, test_u, test_v, simVars, op);


			output_filename = write_file(test_h, "reference_h", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file(test_u, "reference_u", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file(test_v, "reference_v", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;



			output_filename = write_file((prog_h-test_h), "ref_diff_h", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file((prog_u-test_u), "ref_diff_u", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file((prog_v-test_v), "ref_diff_v", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;
		}
	}



	void timestep_do_output()
	{
		write_file_output();

#if 0
		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();

			// Print header
			if (simVars.timecontrol.current_timestep_nr == 0)
			{
//				std::cout << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";
//				std::cout << std::endl;
			}

			std::cout << std::setprecision(simVars.misc.output_floating_point_precision);
			std::cout << std::endl;

			//Print simulation time, energy and pot enstrophy
			std::cout << "DIAGNOSTIC - time, mass, energy, potential_enstrophy ";
			std::cout << simVars.timecontrol.current_simulation_time << "\t";
			std::cout << simVars.diag.total_mass << "\t";
			std::cout << simVars.diag.total_energy << "\t";
			std::cout << simVars.diag.total_potential_enstrophy << std::endl;
		}
#endif

		if (simVars.misc.verbosity > 0)
		{
			std::cout << "prog_h min/max:\t" << prog_h.physical_reduce_min() << ", " << prog_h.physical_reduce_max() << std::endl;
		}


		if (	param_compute_error &&
				simVars.misc.use_nonlinear_equations == 0 &&
				simVars.setup.benchmark_scenario_id == 4
		)
		{
			SphereData test_h(sphereDataConfig);
			SphereData test_u(sphereDataConfig);
			SphereData test_v(sphereDataConfig);

			BenchmarkCombined::setupInitialConditions(test_h, test_u, test_v, simVars, op);

			std::cout << "ERRORS - time, RMS(h,u,v), MAXABS(h,u,v):\t";
			std::cout << simVars.timecontrol.current_simulation_time << "\t";

			std::cout << (test_h-prog_h).physical_reduce_rms() << "\t";
			std::cout << (test_u-prog_u).physical_reduce_rms() << "\t";
			std::cout << (test_v-prog_v).physical_reduce_rms() << "\t";

			std::cout << (test_h-prog_h).physical_reduce_max_abs() << "\t";
			std::cout << (test_u-prog_u).physical_reduce_max_abs() << "\t";
			std::cout << (test_v-prog_v).physical_reduce_max_abs() << std::endl;
		}


		if (simVars.misc.output_each_sim_seconds > 0)
			while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;
	}

public:
	bool timestep_check_output()
	{
		std::cout << "." << std::flush;

		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
			return false;

		timestep_do_output();

		return true;
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


	bool instability_detected()
	{
		double max_abs_value = std::abs(simVars.sim.h0)*2.0+1.0;
		if (
				prog_h.physical_reduce_max_abs() > max_abs_value &&
				simVars.setup.benchmark_scenario_id != 4
		)
		{
			std::cerr << "Instability detected (max abs value of h > " << max_abs_value << ")" << std::endl;
			return true;
		}

		if (prog_h.physical_isAnyNaNorInf())
		{
			std::cerr << "Inf value detected" << std::endl;
			return true;
		}

		return false;
	}



	void run_timestep()
	{
#if SWEET_GUI
		if (simVars.misc.gui_enabled)
			timestep_check_output();
#endif

		// output of time step size
		double o_dt;

		if (simVars.rexi.use_rexi == 0)
		{
			timestepping_explicit.run_rk_timestep(
					this,
					&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					prog_h, prog_u, prog_v,
					o_dt,
					simVars.timecontrol.current_timestep_size,
					simVars.disc.timestepping_runge_kutta_order,
					simVars.timecontrol.current_simulation_time,
					simVars.timecontrol.max_simulation_time
				);
		}
		else if (simVars.rexi.use_rexi == -1)
		{
			SphereData o_prog_h(sphereDataConfig);
			SphereData o_prog_u(sphereDataConfig);
			SphereData o_prog_v(sphereDataConfig);

			timestepping_implicit_swe.solve(
					prog_h*simVars.sim.gravitation, prog_u, prog_v,
					o_prog_h, o_prog_u, o_prog_v
				);

			prog_h = o_prog_h / simVars.sim.gravitation;
			prog_u = o_prog_u;
			prog_v = o_prog_v;
		}
		else
		{
			o_dt = simVars.timecontrol.current_timestep_size;
			assert(o_dt > 0);

			// padding to max simulation time if exceeding the maximum
			if (simVars.timecontrol.max_simulation_time >= 0)
				if (o_dt + simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
					o_dt = simVars.timecontrol.max_simulation_time-simVars.timecontrol.current_simulation_time;

			swe_sphere_rexi.run_timestep_rexi(
					prog_h,
					prog_u,
					prog_v,
					o_dt,
					simVars
				);

			/*
			 * Add implicit viscosity
			 */
			if (simVars.sim.viscosity != 0)
			{
				double scalar = simVars.sim.viscosity*simVars.timecontrol.current_timestep_size;
				double r = simVars.sim.earth_radius;

				/*
				 * (1-dt*visc*D2)p(t+dt) = p(t)
				 */
				prog_h = prog_h.spectral_solve_helmholtz(1.0, -scalar, r);
				prog_u = prog_u.spectral_solve_helmholtz(1.0, -scalar, r);
				prog_v = prog_v.spectral_solve_helmholtz(1.0, -scalar, r);
			}
		}

		// advance time step and provide information to parameters
		simVars.timecontrol.current_timestep_size = o_dt;
		simVars.timecontrol.current_simulation_time += o_dt;
		simVars.timecontrol.current_timestep_nr++;
	}


	// Main routine for method to be used in case of finite differences
	void p_run_euler_timestep_update(
			const SphereData &i_h,	///< prognostic variables
			const SphereData &i_u,	///< prognostic variables
			const SphereData &i_v,	///< prognostic variables

			SphereData &o_h_t,	///< time updates
			SphereData &o_u_t,	///< time updates
			SphereData &o_v_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;

		if (!simVars.misc.use_nonlinear_equations)
		{
			if (!simVars.misc.sphere_use_robert_functions)
			{
				// linear equations
				o_h_t = -(op.div_lon(i_u)+op.div_lat(i_v))*(simVars.sim.h0/simVars.sim.earth_radius);

				o_u_t = -op.grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					o_u_t += f(i_v);
					o_v_t -= f(i_u);
				}
			}
			else
			{
				// use Robert functions for velocity
				// linear equations
				o_h_t = -(op.robert_div_lon(i_u)+op.robert_div_lat(i_v))*(simVars.sim.h0/simVars.sim.earth_radius);

				o_u_t = -op.robert_grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.robert_grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					o_u_t += f(i_v);
					o_v_t -= f(i_u);
				}
			}
		}
		else
		{
			assert(simVars.sim.earth_radius > 0);
			assert(simVars.sim.gravitation);

			if (simVars.misc.sphere_use_robert_functions)
			{
				FatalError("Only non-robert formulation is supported so far for non-linear SWE on sphere!");
				// TODO: rewrite for robert functions
				// TODO: Also initialize velocities correctly
			}

			/*
			 * Height
			 */
			// non-linear equations
			o_h_t = -(op.div_lon(i_h*i_u)+op.div_lat(i_h*i_v))*(1.0/simVars.sim.earth_radius);

			/*
			 * Velocity
			 */
			// linear terms
			o_u_t = -op.grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
			o_v_t = -op.grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

			if (simVars.sim.coriolis_omega != 0)
			{
				o_u_t += f(i_v);
				o_v_t -= f(i_u);
			}

			// non-linear terms
			o_u_t -= (i_u*op.grad_lon(i_u) + i_v*op.grad_lat(i_u))*(1.0/simVars.sim.earth_radius);
			o_v_t -= (i_u*op.grad_lon(i_v) + i_v*op.grad_lat(i_v))*(1.0/simVars.sim.earth_radius);
		}

		assert(simVars.sim.viscosity_order == 2);
		if (simVars.sim.viscosity != 0)
		{
			double scalar = simVars.sim.viscosity/(simVars.sim.earth_radius*simVars.sim.earth_radius);

			o_h_t += op.laplace(i_h)*scalar;
			o_u_t += op.laplace(i_u)*scalar;
			o_v_t += op.laplace(i_v)*scalar;
		}
	}

#if SWEET_GUI


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



	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data
	)
	{
		// request rendering of sphere
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		int id = simVars.misc.vis_id % 4;
		switch (id)
		{
		default:
		case 0:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_h, planeDataConfig);
			break;

		case 1:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(
					prog_u
//					simVars.misc.sphere_use_robert_functions ? prog_u.robert_convertToNonRobert() : prog_u
					, planeDataConfig);
			break;

		case 2:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(
					prog_v
//					simVars.misc.sphere_use_robert_functions ? prog_v.robert_convertToNonRobert() : prog_v
					, planeDataConfig);
			break;

		case 3:
			if (simVars.misc.sphere_use_robert_functions)
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(op.vort(prog_u, prog_v), planeDataConfig);
			else
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(op.robert_vort(prog_u, prog_v), planeDataConfig);
			break;

		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostic values if required
//		update_diagnostics();

		const char* description = "";

		int id = simVars.misc.vis_id % 4;
		switch (id)
		{
		default:
		case 0:
			description = "H";
			break;

		case 1:
			description = "U";
			break;

		case 2:
			description = "V";
			break;

		case 3:
			description = "eta";
			break;
		}

		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string, "Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy,
				viz_plane_data.reduce_max(),
				viz_plane_data.reduce_min()
		);

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

		case 'b':
			render_primitive_id = (render_primitive_id + 1) % 2;
			break;

		case 'c':
			write_file_output();
			break;

#if 0
case 'C':
			// dump data arrays to VTK
			prog_h.file_physical_saveData_vtk("swe_rexi_dump_h.vtk", "Height");
			prog_u.file_physical_saveData_vtk("swe_rexi_dump_u.vtk", "U-Velocity");
			prog_v.file_physical_saveData_vtk("swe_rexi_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_h.file_physical_loadData("swe_rexi_dump_h.csv", simVars.setup.input_data_binary);
			prog_u.file_physical_loadData("swe_rexi_dump_u.csv", simVars.setup.input_data_binary);
			prog_v.file_physical_loadData("swe_rexi_dump_v.csv", simVars.setup.input_data_binary);
			break;
#endif
		}
	}
#endif
};


int main(int i_argc, char *i_argv[])
{

#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

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

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"rexi-use-coriolis-formulation",
			"compute-error",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 1;
	simVars.bogus.var[1] = 1;
	simVars.bogus.var[2] = 1;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		std::cout << "	--compute-error [0/1]	Output errors (if available, default: 1)" << std::endl;
		std::cout << "	--rexi-use-coriolis-formulation [0/1]	Use Coriolisincluding  solver for REXI (default: 1)" << std::endl;
		std::cout << "	--benchmark-geostrophic-balance-freq-multiplier [int]	(default: 1)" << std::endl;
		return -1;
	}

	param_rexi_use_coriolis_formulation = simVars.bogus.var[0];
	assert (param_rexi_use_coriolis_formulation == 0 || param_rexi_use_coriolis_formulation == 1);

	param_compute_error = simVars.bogus.var[1];
	param_geostr_balance_freq_multiplier = simVars.bogus.var[2];

	sphereDataConfigInstance.setupAutoPhysicalSpace(
					simVars.disc.res_spectral[0],
					simVars.disc.res_spectral[1],
					&simVars.disc.res_physical[0],
					&simVars.disc.res_physical[1]
			);

	sphereDataConfigInstanceExt.setupAdditionalModes(
			&sphereDataConfigInstance,
			simVars.rexi.rexi_use_extended_modes,
			simVars.rexi.rexi_use_extended_modes
		);

#if SWEET_GUI
	planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical);
#endif

	std::ostringstream buf;
	buf << std::setprecision(14);

	if (simVars.setup.benchmark_scenario_id == 1)
	{
		/// Setup Galewski parameters
		simVars.sim.coriolis_omega = 7.292e-5;
		simVars.sim.gravitation = 9.80616;
		simVars.sim.earth_radius = 6.37122e6;
		simVars.sim.h0 = 10000.0;

		simVars.misc.output_time_scale = 1.0/(60.0*60.0);
	}


#if SWEET_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// only start simulation and time stepping for first rank
	if (rank == 0)
#endif
	{
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

#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			VisSweet<SimulationInstance> visSweet(simulationSWE);
			delete simulationSWE;
		}
		else
#endif
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			//Setting initial conditions and workspace - in case there is no GUI

			// reset already triggered
			//simulationSWE->reset();

			//Time counter
			Stopwatch time;

			//Diagnostic measures at initial stage
			//double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

			// Initialize diagnostics
			if (simVars.misc.verbosity > 0)
			{
//				simulationSWE->update_diagnostics();
#if 0
				diagnostics_energy_start = simVars.diag.total_energy;
				diagnostics_mass_start = simVars.diag.total_mass;
				diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
#endif
			}

#if SWEET_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			//Start counting time
			time.reset();


			// Main time loop
			while(true)
			{
				simulationSWE->timestep_check_output();

				// Stop simulation if requested
				if (simulationSWE->should_quit())
					break;

				// Main call for timestep run
				simulationSWE->run_timestep();

				// Instability
				if (simulationSWE->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}
			}

			// Always write or overwrite final time step output!
			simulationSWE->timestep_do_output();

			// Stop counting time
			time.stop();

			double seconds = time();

			// End of run output results
			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;

			delete simulationSWE;
		}
	}
#if SWEET_MPI
	else
	{
		if (simVars.rexi.use_rexi == 1)
		{
			SWE_Sphere_REXI rexiSWE;

			/*
			 * Setup our little dog REXI
			 */
			rexiSWE.setup(
					simVars.rexi.rexi_h,
					simVars.rexi.rexi_M,
					simVars.rexi.rexi_L,

					sphereDataConfig,
					&simVars.sim,
					simVars.timecontrol.current_timestep_size,

					simVars.rexi.rexi_use_half_poles,
					simVars.misc.sphere_use_robert_functions,
					simVars.rexi.rexi_use_extended_modes,
					simVars.rexi.rexi_normalization,
					param_rexi_use_coriolis_formulation
				);

			bool run = true;

			SphereData prog_h(sphereDataConfig);
			SphereData prog_u(sphereDataConfig);
			SphereData prog_v(sphereDataConfig);

			// initialize with dummy data
			prog_h.spectral_set_zero();
			prog_u.spectral_set_zero();
			prog_v.spectral_set_zero();

			MPI_Barrier(MPI_COMM_WORLD);

			while (run)
			{
				// REXI time stepping
				run = rexiSWE.run_timestep_rexi(
						prog_h, prog_u, prog_v,
						-simVars.sim.CFL,
						simVars
				);
			}
		}
	}
#endif


#if SWEET_MPI
	if (simVars.rexi.use_rexi == 1)
	{
		// synchronize REXI
		if (rank == 0)
			SWE_Sphere_REXI::MPI_quitWorkers(sphereDataConfig);
	}

	MPI_Finalize();
#endif

	return 0;
}
