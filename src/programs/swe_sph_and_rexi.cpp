/*
 * AppTestSWE.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_TESTSWE_HPP_
#define SRC_TESTSWE_HPP_

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
	#include <sweet/plane/PlaneDataConfig.hpp>
	#include <sweet/plane/PlaneData.hpp>
	#include <sweet/Convert_SphereData_To_PlaneData.hpp>
#endif

#include <benchmarks_sphere/BenchmarkGalewsky.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataTimesteppingRK.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphBandedMatrix.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>


#include "swe_sphere_rexi/SWE_Sphere_REXI.hpp"



SimulationVariables simVars;

//Diagnostic measures at initial stage
//double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

// Plane data config
SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;


#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif

/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */
bool param_rexi_use_coriolis_formulation = false;



class SimulationInstance
{
public:
	SphereOperators op;
	SphereOperatorsComplex opComplex;

	// Runge-Kutta stuff
	SphereDataTimesteppingRK timestepping;

	SWE_Sphere_REXI swe_sphere_rexi;


	SphereData prog_h;
	SphereData prog_u;
	SphereData prog_v;


	BenchmarkGalewsky benchmarkGalewsky;

	REXI rexi;

#if SWEET_GUI
	PlaneData viz_plane_data;
#endif


public:
	SimulationInstance()	:
		prog_h(sphereDataConfig),
		prog_u(sphereDataConfig),
		prog_v(sphereDataConfig),
		benchmarkGalewsky(::simVars)

#if SWEET_GUI
		,viz_plane_data(planeDataConfig)
#endif
	{
		reset();
	}



	void write_file_output()
	{
		char buffer[1024];

		std::cout << "Simulation time: " << simVars.timecontrol.current_simulation_time << std::endl;

		sprintf(buffer, "prog_h_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (simVars.setup.benchmark_scenario_id == 0)
			prog_h.file_physical_writeFile_lon_pi_shifted(buffer);
		else
			prog_h.physical_write_file(buffer);
		std::cout << buffer << " (min: " << prog_h.reduce_min() << ", max: " << prog_h.reduce_max() << ")" << std::endl;

		sprintf(buffer, "prog_u_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (simVars.setup.benchmark_scenario_id == 0)
			prog_u.file_physical_writeFile_lon_pi_shifted(buffer);
		else
			prog_u.physical_write_file(buffer);
		std::cout << buffer << std::endl;

		sprintf(buffer, "prog_v_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (simVars.setup.benchmark_scenario_id == 0)
			prog_v.file_physical_writeFile_lon_pi_shifted(buffer);
		else
			prog_v.physical_write_file(buffer);
		std::cout << buffer << std::endl;

		sprintf(buffer, "prog_eta_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		SphereData vort = op.vort(prog_u, prog_v)/simVars.sim.earth_radius;
		if (simVars.setup.benchmark_scenario_id == 0)
			vort.file_physical_writeFile_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			vort.physical_write_file(buffer);
		std::cout << buffer << std::endl;
	}



	void setup_initial_conditions_gaussian(
			double i_center_lat = M_PI/3,
			double i_center_lon = M_PI/3
	)
	{
		double exp_fac = 10.0;

		double center_lat = i_center_lat;
		double center_lon = i_center_lon;

		auto initial_condition_h = [&](double lon, double mu, double &o_data)
		{
			// https://en.wikipedia.org/wiki/Great-circle_distance
			// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
			// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			o_data = exp(-d*d*exp_fac)*0.1*simVars.sim.h0 + simVars.sim.h0;
		};

		prog_h.physical_update_lambda_gaussian_grid(initial_condition_h);
		prog_u.physical_set_zero();
		prog_v.physical_set_zero();
	}



	SphereData f(SphereData i_sphData)
	{
		return op.mu(i_sphData*2.0*simVars.sim.coriolis_omega);
	}


	void update_diagnostics()
	{
	}

	void reset()
	{
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


		if (simVars.setup.benchmark_scenario_id == 0)
		{
			setup_initial_conditions_gaussian(M_PI/3.0, M_PI/3.0);
		}
		else if (simVars.setup.benchmark_scenario_id == 1)
		{
			/// setup galewski parameters
			simVars.sim.coriolis_omega = 7.292e-5;
			simVars.sim.gravitation = 9.80616;
			simVars.sim.earth_radius = 6.37122e6;
			simVars.sim.h0 = 10000.0;

			benchmarkGalewsky.setup_initial_h(prog_h);
//			prog_h.spat_set_zero();
			benchmarkGalewsky.setup_initial_h_add_bump(prog_h);

			benchmarkGalewsky.setup_initial_u(prog_u);
			benchmarkGalewsky.setup_initial_v(prog_v);
		}
		else if (simVars.setup.benchmark_scenario_id == 2 || simVars.setup.benchmark_scenario_id == 3)
		{
#if 0
			if (simVars.timecontrol.current_timestep_size <= 0)
			{
				std::cout << "Timestep size not positive" << std::endl;
				assert(false);
				exit(1);
			}
#endif

			if (simVars.setup.benchmark_scenario_id == 2)
			{
				setup_initial_conditions_gaussian(0);
			}
			else if (simVars.setup.benchmark_scenario_id == 3)
			{
				setup_initial_conditions_gaussian(M_PI/3.0);
				//setup_initial_conditions_gaussian(M_PI, M_PI);
//				setup_initial_conditions_gaussian(M_PI*0.5, M_PI*0.5);
//				setup_initial_conditions_gaussian(-M_PI/3.0);
			}
		}



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
		std::cout << " + Use REXI: " << simVars.rexi.use_rexi << std::endl;
		std::cout << " + REXI h: " << simVars.rexi.rexi_h << std::endl;
		std::cout << " + REXI M: " << simVars.rexi.rexi_M << std::endl;
		std::cout << " + REXI use half poles: " << simVars.rexi.rexi_use_half_poles << std::endl;
		std::cout << " + REXI additional modes: " << simVars.rexi.rexi_use_extended_modes << std::endl;
		std::cout << std::endl;
		std::cout << " + RK order: " << simVars.disc.timestepping_runge_kutta_order << std::endl;
		std::cout << " + timestep size: " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << " + output timestep size: " << simVars.misc.output_each_sim_seconds << std::endl;

		std::cout << std::endl;

		if (simVars.rexi.use_rexi)
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
					param_rexi_use_coriolis_formulation
				);
		}
	}


public:
	bool timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		std::cout << "." << std::flush;

		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
			return false;

		write_file_output();

		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();

			// Print header
			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "T\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";

				//if ((simVars.setup.scenario >= 0 && simVars.setup.scenario <= 4) || simVars.setup.scenario == 13)
				o_ostream << "\tDIFF_H0\tDIFF_U0\tDIFF_V0";

				if (simVars.misc.use_nonlinear_equations==0){
					o_ostream << "\tANAL_DIFF_RMS_P\tANAL_DIFF_RMS_U\tANAL_DIFF_RMS_V";
					o_ostream << "\tANAL_DIFF_MAX_P\tANAL_DIFF_MAX_U\tANAL_DIFF_MAX_V";
				}
				o_ostream << std::endl;
			}

			//Print simulation time, energy and pot enstrophy
			o_ostream << std::setprecision(8) << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_mass << "\t" << simVars.diag.total_energy << "\t" << simVars.diag.total_potential_enstrophy;
		}


		if (simVars.misc.output_each_sim_seconds > 0)
			while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;

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
		if (prog_h.reduce_abs_max() > max_abs_value)
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
		if (simVars.rexi.use_rexi == false)
		{
			double o_dt;
			timestepping.run_rk_timestep(
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
		else
		{
			double o_dt = simVars.timecontrol.current_timestep_size;
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

			// advance time step and provide information to parameters
			simVars.timecontrol.current_timestep_size = o_dt;
			simVars.timecontrol.current_simulation_time += o_dt;
			simVars.timecontrol.current_timestep_nr++;
		}

#if SWEET_GUI
		timestep_output();
#endif
	}

#if 0
	void run()
	{
		reset();


//		// This class is only used in case of added modes
//		SphereDataConfig sphereDataConfigRexiAddedModes;

		// Pointer to SPH configuration for REXI computations
//		SphereDataConfig *sphereDataConfigRexi = nullptr;


		// Main time loop
		while(true)
		{

		}

		if (simVars.rexi.use_rexi == false)
		{
			simVars.timecontrol.current_simulation_time = 0;
			while (!should_quit())
			{
				if (simVars.timecontrol.current_simulation_time >= simVars.misc.output_next_sim_seconds)
				{
					std::cout << std::endl;
					write_file_output();

					simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;
					if (simVars.misc.output_next_sim_seconds < simVars.timecontrol.current_simulation_time)
						simVars.misc.output_next_sim_seconds = simVars.timecontrol.current_simulation_time;
				}

				double o_dt;
				timestepping.run_rk_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
						prog_h, prog_u, prog_v,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_runge_kutta_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);

				std::cout << "." << std::flush;

#if 1
				double max_abs_value = std::abs(simVars.sim.h0)*2.0+1.0;
				if (prog_h.reduce_abs_max() > max_abs_value)
				{
					std::cerr << "Instability detected (max abs value of h > " << max_abs_value << ")" << std::endl;
					assert(false);
					exit(1);
				}
#endif
				if (prog_h.physical_isAnyNaNorInf())
				{
					std::cerr << "Instability detected (NaN value in H)" << std::endl;
					assert(false);
					exit(1);
				}

				if (simVars.timecontrol.current_simulation_time >= simVars.timecontrol.max_simulation_time)
					break;

				advance_simulation_time();
			}

			std::cout << std::endl;
			write_file_output();
			std::cout << std::endl;
		}
		else
		{
#if 1
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
					param_rexi_use_coriolis_formulation
				);

#else
			if (simVars.rexi.rexi_use_extended_modes == 0)
			{
				sphereDataConfigRexi = sphereDataConfig;
			}
			else
			{
				// Add modes only along latitude since these are the "problematic" modes
				sphereDataConfigRexiAddedModes.setupAdditionalModes(
						sphereDataConfig,
						simVars.rexi.rexi_use_extended_modes,	// TODO: Extend SPH wrapper to also support m != n to set this guy to 0
						simVars.rexi.rexi_use_extended_modes
				);
				sphereDataConfigRexi = &sphereDataConfigRexiAddedModes;
			}

			rexi.setup(simVars.rexi.rexi_h, simVars.rexi.rexi_M, 0, simVars.rexi.rexi_use_half_poles);

			std::cout << "REXI poles: " << rexi.alpha.size() << std::endl;

			SphereData tmp_prog_phi(sphereDataConfigRexi);
			SphereData tmp_prog_u(sphereDataConfigRexi);
			SphereData tmp_prog_v(sphereDataConfigRexi);

			SphereData accum_prog_phi(sphereDataConfigRexi);
			SphereData accum_prog_u(sphereDataConfigRexi);
			SphereData accum_prog_v(sphereDataConfigRexi);

			std::vector<SWERexi_SPHRobert> rexiSPHRobert_vector;
			std::vector<SWERexi_SPH> rexiSPH_vector;

			bool use_rexi_preallocaation = true;

			if (use_rexi_preallocaation)
			{
				if (simVars.misc.sphere_use_robert_functions)
					rexiSPHRobert_vector.resize(rexi.alpha.size());
				else
					rexiSPH_vector.resize(rexi.alpha.size());


				for (std::size_t i = 0; i < rexi.alpha.size(); i++)
				{
					std::complex<double> &alpha = rexi.alpha[i];
					std::complex<double> &beta = rexi.beta_re[i];

					if (simVars.misc.sphere_use_robert_functions)
					{
						rexiSPHRobert_vector[i].setup(
								sphereDataConfigRexi,
								alpha,
								beta,
								simVars.sim.earth_radius,
								simVars.sim.coriolis_omega, //*inv_sqrt_avg_geopo,
								simVars.sim.h0*simVars.sim.gravitation,
								simVars.timecontrol.current_timestep_size, //*sqrt_avg_geopo
								param_rexi_use_coriolis_formulation
						);
					}
					else
					{
						rexiSPH_vector[i].setup(
								sphereDataConfigRexi,
								alpha,
								beta,
								simVars.sim.earth_radius,
								simVars.sim.coriolis_omega, //*inv_sqrt_avg_geopo,
								simVars.sim.h0*simVars.sim.gravitation,
								simVars.timecontrol.current_timestep_size, //*sqrt_avg_geopo
								param_rexi_use_coriolis_formulation
						);
					}
				}
			}
#endif
			simVars.timecontrol.current_simulation_time = 0;
			while (simVars.timecontrol.current_simulation_time <= simVars.timecontrol.max_simulation_time)
			{
				if (simVars.timecontrol.current_simulation_time >= simVars.misc.output_next_sim_seconds)
				{
					write_file_output();

					simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;
					if (simVars.misc.output_next_sim_seconds < simVars.timecontrol.current_simulation_time)
						simVars.misc.output_next_sim_seconds = simVars.timecontrol.current_simulation_time;
				}


#if 1

				swe_sphere_rexi.run_timestep_rexi(prog_h, prog_u, prog_v, simVars.timecontrol.current_timestep_size, simVars);

#else

				{
					// convert to geopotential
					SphereData prog_phi_rexi(sphereDataConfigRexi);
					SphereData prog_u_rexi(sphereDataConfigRexi);
					SphereData prog_v_rexi(sphereDataConfigRexi);

					if (simVars.rexi.rexi_use_extended_modes == 0)
					{
						prog_phi_rexi = prog_h*simVars.sim.gravitation;
						prog_u_rexi = prog_u;
						prog_v_rexi = prog_v;
					}
					else
					{
						(prog_h*simVars.sim.gravitation).spectral_copyToDifferentModes(prog_phi_rexi);
						prog_u.spectral_copyToDifferentModes(prog_u_rexi);
						prog_v.spectral_copyToDifferentModes(prog_v_rexi);
					}

					accum_prog_phi.physical_set_zero();
					accum_prog_u.physical_set_zero();
					accum_prog_v.physical_set_zero();

					for (std::size_t i = 0; i < rexi.alpha.size(); i++)
					{
						std::complex<double> &alpha = rexi.alpha[i];
						std::complex<double> &beta = rexi.beta_re[i];

						if (simVars.misc.sphere_use_robert_functions)
						{
							if (use_rexi_preallocaation)
							{
								rexiSPHRobert_vector[i].solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
							else
							{
								SWERexi_SPHRobert rexiSPHRobert;
								rexiSPHRobert.setup(
										sphereDataConfig,
										alpha,
										beta,
										simVars.sim.earth_radius,
										simVars.sim.coriolis_omega, //*inv_sqrt_avg_geopo,
										simVars.sim.h0*simVars.sim.gravitation,
										simVars.timecontrol.current_timestep_size, //*sqrt_avg_geopo
										param_rexi_use_coriolis_formulation
								);

								rexiSPHRobert.solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
						}
						else
						{
							if (use_rexi_preallocaation)
							{
								rexiSPH_vector[i].solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
							else
							{
								SWERexi_SPH rexiSPH;
								rexiSPH.setup(
										sphereDataConfig,
										alpha,
										beta,
										simVars.sim.earth_radius,
										simVars.sim.coriolis_omega,
										simVars.sim.h0*simVars.sim.gravitation,
										simVars.timecontrol.current_timestep_size,
										param_rexi_use_coriolis_formulation
								);

								rexiSPH.solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
						}

						accum_prog_phi += tmp_prog_phi;
						accum_prog_u += tmp_prog_u;
						accum_prog_v += tmp_prog_v;
					}

					if (simVars.rexi.rexi_use_extended_modes == 0)
					{
						prog_h = accum_prog_phi*(1.0/simVars.sim.gravitation);
						prog_u = accum_prog_u;
						prog_v = accum_prog_v;
					}
					else
					{
						(accum_prog_phi*(1.0/simVars.sim.gravitation)).spectral_copyToDifferentModes(prog_h);
						accum_prog_u.spectral_copyToDifferentModes(prog_u);
						accum_prog_v.spectral_copyToDifferentModes(prog_v);
					}
				}

#endif

				std::cout << "." << std::flush;

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

#if 1
				double max_abs_value = std::abs(simVars.sim.h0)*2.0+1.0;
				if (prog_h.reduce_abs_max() > max_abs_value)
				{
					std::cerr << "Instability detected (max abs value of h > " << max_abs_value << ")" << std::endl;
					assert(false);
					exit(1);
				}
#endif
				if (simVars.timecontrol.current_simulation_time >= simVars.timecontrol.max_simulation_time)
					break;

				advance_simulation_time();
			}

			write_file_output();
			std::cout << std::endl;


			if (prog_h.physical_isAnyNaNorInf())
			{
				std::cerr << "Instability detected (NaN value in H)" << std::endl;
				assert(false);
				exit(1);
			}
		}
	}
#endif


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
			double *o_aspect_ratio
	)
	{
		int id = simVars.misc.vis_id % 4;
		switch (id)
		{
		default:
		case 0:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_h, planeDataConfig);
			break;

		case 1:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_u, planeDataConfig);
			break;

		case 2:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_v, planeDataConfig);
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
		update_diagnostics();

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
		sprintf(title_string, "Time (days): %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
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
#if 0
		case 'c':
			// dump data arrays
			prog_h.file_physical_saveData_ascii("swe_rexi_dump_h.csv");
			prog_u.file_physical_saveData_ascii("swe_rexi_dump_u.csv");
			prog_v.file_physical_saveData_ascii("swe_rexi_dump_v.csv");
			break;

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
	MemBlockAlloc::setup();

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"rexi-use-coriolis-formulation",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 0;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}

	param_rexi_use_coriolis_formulation = simVars.bogus.var[0];
	assert (param_rexi_use_coriolis_formulation == 0 || param_rexi_use_coriolis_formulation == 1);


	sphereDataConfigInstance.setupAutoPhysicalSpace(
					simVars.disc.res_spectral[0],
					simVars.disc.res_spectral[1],
					&simVars.disc.res_physical[0],
					&simVars.disc.res_physical[1]
			);


#if SWEET_GUI
	planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical);
#endif


	std::ostringstream buf;
	buf << std::setprecision(14);

#if 0
	SimulationInstance test_swe(sphereDataConfig);
	test_swe.run();
#else

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

			simulationSWE->reset();

			//Time counter
			Stopwatch time;

			//Diagnostic measures at initial stage
			//double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

			// Initialize diagnostics
			if (simVars.misc.verbosity > 0)
			{
				simulationSWE->update_diagnostics();
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
				if (simulationSWE->timestep_output(buf))
				{
					// string output data

					std::string output = buf.str();
					buf.str("");

					// This is an output printed on screen or buffered to files if > used
					std::cout << output;
				}

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

			// Stop counting time
			time.stop();

			double seconds = time();

			// End of run output results
			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;
#if 0
			if (param_timestepping_mode != 0)
				std::cout << "REXI alpha.size(): " << simulationSWE->swe_plane_rexi.rexi.alpha.size() << std::endl;

			if (simVars.misc.verbosity > 0)
			{
				std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((simVars.diag.total_energy-diagnostics_energy_start)/diagnostics_energy_start) << std::endl;
				std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
				std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((simVars.diag.total_potential_enstrophy-diagnostics_potential_entrophy_start)/diagnostics_potential_entrophy_start) << std::endl;

				if (param_compute_error)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark_diff_h << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
				}
			}

			if (param_compute_error && param_nonlinear==0)
			{
				simulationSWE->compute_errors();
				std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << simulationSWE->benchmark_analytical_error_rms_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << simulationSWE->benchmark_analytical_error_rms_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << simulationSWE->benchmark_analytical_error_rms_v << std::endl;

				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << simulationSWE->benchmark_analytical_error_maxabs_h << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << simulationSWE->benchmark_analytical_error_maxabs_u << std::endl;
				std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << simulationSWE->benchmark_analytical_error_maxabs_v << std::endl;
			}

#endif

			delete simulationSWE;
		}
	}
#if SWEET_MPI
	else
	{
		if (param_timestepping_mode == 1)
		{
			SWE_Plane_REXI rexiSWE;

			/*
			 * Setup our little dog REXI
			 */
			rexiSWE.setup(
					simVars.rexi.rexi_h,
					simVars.rexi.rexi_m,
					simVars.rexi.rexi_l,
					simVars.disc.res_physical,
					simVars.sim.domain_size,
					simVars.rexi.rexi_half,
					simVars.rexi.rexi_use_spectral_differences_for_complex_array,
					simVars.rexi.rexi_helmholtz_solver_id,
					simVars.rexi.rexi_helmholtz_solver_eps
				);

			bool run = true;

			PlaneData prog_h(planeDataConfig);
			PlaneData prog_u(planeDataConfig);
			PlaneData prog_v(planeDataConfig);

			PlaneOperators op(simVars.disc.res_physical, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

			MPI_Barrier(MPI_COMM_WORLD);

			while (run)
			{
				// REXI time stepping
				run = rexiSWE.run_timestep_rexi(
						prog_h, prog_u, prog_v,
						-simVars.sim.CFL,
						op,
						simVars
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
			SWE_Plane_REXI::MPI_quitWorkers(planeDataConfig);
	}

	MPI_Finalize();
#endif

#endif
	return 0;
}



#endif /* SRC_TESTSWE_HPP_ */
