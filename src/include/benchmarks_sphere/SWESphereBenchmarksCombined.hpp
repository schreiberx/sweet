/*
 * BenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *	  Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SPHERE_SWESPHEREBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_SPHERE_SWESPHEREBENCHMARKSCOMBINED_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>
//#include <benchmarks_sphere/BenchmarkGalewsky.hpp>
#include <benchmarks_sphere/BenchmarkGaussianDam.hpp>
#include <benchmarks_sphere/BenchmarkFlowOverMountain.hpp>

#if SWEET_MPI
	#include <mpi.h>
#endif

class SWESphereBenchmarksCombined
{
public:
	// Simulation variables
	SimulationVariables *simVars;

	// Sphere operators
	SphereOperators *op;

	// plane or sphere data config
	const void* ext_forces_data_config;


	SWESphereBenchmarksCombined()	:
		simVars(nullptr),
		op(nullptr),
		ext_forces_data_config(nullptr)
	{
	}


public:
	static void printAvailableBenchmarks()
	{
		std::cout << std::endl;
		std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
		std::cout << "Available benchmark scenarios (DEPRECATED):" << std::endl;
		std::cout << "	0: Dummy" << std::endl;
		std::cout << "	2: Use Gaussian bump initial conditions (0, pi/3, exp=10)" << std::endl;
		std::cout << "	3: Use Gaussian bump initial conditions (pi/3, pi/3, exp=10)" << std::endl;
		std::cout << "	4: Use Gaussian bump initial conditions (-pi/2, pi/4, exp=10)" << std::endl;
		std::cout << "	4: Use Gaussian bump initial conditions (-pi/2, pi/4, exp=10)" << std::endl;
		std::cout << "	5: Use Gaussian bump initial conditions (-pi/2, pi/4, exp=10)" << std::endl;
		std::cout << "	6: Use Gaussian bump initial conditions (pi/3, pi/3, exp=20)" << std::endl;
		std::cout << "	7: Use Gaussian bump initial conditions (-pi, pi/4, exp=100)" << std::endl;
		std::cout << "	9: Combination of Gaussian bumps" << std::endl;

		std::cout << "	10: Williamson Geostrophic balance test case, unit sphere" << std::endl;
		std::cout << "	11: Williamson Geostrophic balance test case, earth scale" << std::endl;

		std::cout << "	21: Williamson Advection test benchmark 1 (DIV formulation)" << std::endl;
		std::cout << "	22: Williamson Advection test benchmark 1 (GRAD formulation)" << std::endl;

		std::cout << "	50: Swartztrauber 2004 Gaussian breaking dam" << std::endl;
		std::cout << "	200: h=h0, u=0, v=0" << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "Benchmark scenario by name (NEW):" << std::endl;
		std::cout << "     'flat': Constant height and zero velocity" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #1:" << std::endl;
		std::cout << "     'adv_cosine_bell': Advection test case of cosine bell" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #1 (alternative):" << std::endl;
		std::cout << "     'adv_gauss_bump': Advection test case of gaussian bump" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #2 [TODO: Check]:" << std::endl;
		std::cout << "     'geostrophic_balance': Geostrophic balance, one wave (standard)" << std::endl;
		std::cout << "     'geostrophic_balance_[N]': Geostrophic balance, with [N] waves" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #3:" << std::endl;
		std::cout << "     [NOT YET IMPLEMENTED]" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #4:" << std::endl;
		std::cout << "     [NOT YET IMPLEMENTED]" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #5:" << std::endl;
		std::cout << "     'williamson5'/'flow_over_mountain': Flow over mountain benchmark" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #6:" << std::endl;
		std::cout << "     'williamson6'/'rossby_haurwitz_wave': Rossby Haurwitz wave" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #7:" << std::endl;
		std::cout << "     [NOT YET IMPLEMENTED]" << std::endl;
		std::cout << std::endl;
		/*
		 * Williamson #3 test case is only simliar to Galewsky, but not identical
		 */
		std::cout << "  'galewsky': Galwesky benchmark" << std::endl;
		std::cout << "  'galewsky_nobump': Galwesky benchmark without any bump" << std::endl;
		std::cout << "  'galewsky_nosetparams': Galwesky benchmark without setting parameters" << std::endl;
		std::cout << std::endl;
	}



	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance(
			SphereData &o_phi,
			SphereData &i_vort,
			SphereData &i_div,
			// Don't use robert formulation per default! Otherwise bugs in here can cancel out bugs in the Robert time stepper
			bool i_use_robert_formulation = false
	)
	{
		/*
		 * Setup Coriolis effect
		 */
		SphereDataPhysical f(o_phi.sphereDataConfig);
		f.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = 2.0*simVars->sim.coriolis_omega*mu;
			}
		);

		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereDataPhysical u(o_phi.sphereDataConfig);
		SphereDataPhysical v(o_phi.sphereDataConfig);

		if (i_use_robert_formulation)
			op->robert_vortdiv_to_uv(i_vort, i_div, u, v);
		else
			op->vortdiv_to_uv(i_vort, i_div, u, v);

		SphereDataPhysical vrtg = i_vort.getSphereDataPhysical();

		SphereDataPhysical tmpg1 = u*(vrtg+f);
		SphereDataPhysical tmpg2 = v*(vrtg+f);

		SphereData tmpspec1(o_phi.sphereDataConfig);
		SphereData tmpspec2(o_phi.sphereDataConfig);

		if (i_use_robert_formulation)
		  op->robert_uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);
		else
		  op->uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		tmpspec2 = 0.5*(u*u+v*v);

		if (i_use_robert_formulation)
			tmpspec2 = tmpspec2.robert_convertToNonRobertSquared();

		SphereData phispec = op->inv_laplace(tmpspec1) - tmpspec2;

		o_phi = phispec + simVars->sim.h0*simVars->sim.gravitation;
	}


public:
	void setupTopography()
	{
		if (simVars == nullptr)
			FatalError("Benchmarks are not yet initialized");

		if (simVars->setup.benchmark_name == "flow_over_mountain")
		{
			// set the topography flag to true
			simVars->sim.use_topography = true;

			// setup the parameters for the flow-over-mountain test case
			const double R			= M_PI/9.;
			const double h_topo_0	 = 2000.;
			const double i_center_lon = 3.*M_PI/2.;
			const double i_center_lat = M_PI/6.;

			// initialize the topography
			simVars->sim.h_topo.physical_set_zero();

			// setup the topography vector
			BenchmarkFlowOverMountain::setup_topography(
					simVars->sim.h_topo,
					*simVars,
					R,
					h_topo_0,
					i_center_lon,
					i_center_lat
			);

			simVars->sim.h_topo.spectral_truncate();
		}
		else
		{
			// set the topography flag to false
			simVars->sim.use_topography = false;
		}
	}



	void setup(
			SimulationVariables &io_simVars,
			SphereOperators &io_op
	)
	{
		simVars = &io_simVars;
		op = &io_op;
	}


	void setupInitialConditions(
			SphereData &o_phi,
			SphereData &o_vort,
			SphereData &o_div
	)
	{
		if (simVars == nullptr)
			FatalError("Benchmarks are not yet initialized");

		simVars->misc.output_time_scale = 1.0/(60.0*60.0);

		if (simVars->setup.benchmark_name != "")
		{
			if (simVars->setup.benchmark_name == "gaussian_bump_advection")
			{
				/*
				 * Advection benchmark with a time-varying velocity field
				 */

				auto callback_external_forces_advection_field =
						[](
								int i_field_id,
								double i_simulation_timestamp,
								void* o_data_void,			/// planedata or spheredata
								void* o_data_user_void		/// user data (pointer to this class)
				)
				{
					SphereData* o_sphere_data = (SphereData*)o_data_void;
					SWESphereBenchmarksCombined* s = (SWESphereBenchmarksCombined*)o_data_user_void;

					if (i_field_id == 1)
					{
						double a = s->simVars->sim.earth_radius;
						double alpha = s->simVars->setup.advection_rotation_angle;
						double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

						double r;
						if (s->simVars->sim.advection_velocity[2] == 0)
							r = 0;
						else
							r = i_simulation_timestamp/s->simVars->sim.advection_velocity[2]*2.0*M_PI;

						// apply rotation
						//alpha += r;

						// change velocity
						u0 = u0*(1.0 + std::cos(r));

						SphereData stream_function(o_sphere_data->sphereDataConfig);

						stream_function.physical_update_lambda(
							[&](double i_lon, double i_lat, double &io_data)
							{
								double i_theta = i_lat;
								double i_lambda = i_lon;

								io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
							}
						);

						*o_sphere_data = s->op->laplace(stream_function);
					}
					else if (i_field_id == 2)
					{
						o_sphere_data->spectral_set_zero();
					}
					else
					{
						FatalError("Non-existing external field requested!");
					}
					return;
				};


				/*
				 * Gaussian bump advection with changing velocity field!
				 */
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					#if SWEET_MPI
						int mpi_rank;
						MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
						if (mpi_rank == 0)
					#endif
					{
						std::cout << "!!! WARNING !!!" << std::endl;
						std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
						std::cout << "!!! WARNING !!!" << std::endl;
					}
				}

				simVars->sim.coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.earth_radius = 6.37122e6;
				simVars->sim.h0 = 1000.0;

				op->setup(o_phi.sphereDataConfig, simVars->sim.earth_radius);

				double lambda_c = 3.0*M_PI/2.0;
				double theta_c = 0.0;

				o_phi.physical_update_lambda(
					[&](double i_lambda, double i_theta, double &io_data)
				{
						double d = std::acos(
								std::sin(theta_c)*std::sin(i_theta) +
								std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
						);

						double i_exp_fac = 20.0;
						io_data = std::exp(-d*d*i_exp_fac)*0.1*simVars->sim.h0;

						io_data *= simVars->sim.gravitation;
					}
				);


				if (simVars->sim.advection_velocity[2] != 0)
				{
					// backup data config
					ext_forces_data_config = o_phi.sphereDataConfig;

					// set callback
					simVars->sim.getExternalForcesCallback = callback_external_forces_advection_field;

					// set user data to this class
					simVars->sim.getExternalForcesUserData = this;
				}

				// setup velocities with initial time stamp
				callback_external_forces_advection_field(1, simVars->timecontrol.current_simulation_time, &o_vort, this);
				callback_external_forces_advection_field(2, simVars->timecontrol.current_simulation_time, &o_div, this);
			}
			else if (
					simVars->setup.benchmark_name == "williamson1"		||
					simVars->setup.benchmark_name == "adv_cosine_bell"
			)
			{
				/*
				 * Advection test case
				 * See Williamson test case, eq. (77), (78), (79)
				 */

				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				simVars->sim.coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.earth_radius = 6.37122e6;
				simVars->sim.h0 = 1000.0;

				// reset operator
				op->setup(o_phi.sphereDataConfig, simVars->sim.earth_radius);

				double lambda_c = 3.0*M_PI/2.0;
				double theta_c = 0.0;
				double a = simVars->sim.earth_radius;

				double R = a/3.0;
				double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

				o_phi.physical_update_lambda(
					[&](double i_lambda, double i_theta, double &io_data)
				{
						double r = a * std::acos(
								std::sin(theta_c)*std::sin(i_theta) +
								std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
						);

						if (r < R)
							io_data = simVars->sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
						else
							io_data = 0;

						io_data *= simVars->sim.gravitation;
					}
				);

				SphereData stream_function(o_phi.sphereDataConfig);

				stream_function.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						double alpha = simVars->setup.advection_rotation_angle;

						io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
					}
				);

				o_vort = op->laplace(stream_function);
				o_div.spectral_set_zero();

	//			std::cout << "advection_rotation_angle: " << simVars->setup.advection_rotation_angle << std::endl;
			}
			else if (
					simVars->setup.benchmark_name == "williamson1b"		||
					simVars->setup.benchmark_name == "adv_gauss_bump"
			)
			{
				/*
				 * Alterative to original Williamson #1 advection test case which is based on a Gaussian bell instead of a cosine bell.
				 * This allows to test for L_inf convergence.
				 */

				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				simVars->sim.coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.earth_radius = 6.37122e6;
				simVars->sim.h0 = 1000.0;

				op->setup(o_phi.sphereDataConfig, simVars->sim.earth_radius);

				double lambda_c = 3.0*M_PI/2.0;
				double theta_c = 0.0;
				double a = simVars->sim.earth_radius;

				//double R = a/3.0;
				double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);
				//double u0 = (2.0*M_PI*a*1000.0);

				o_phi.physical_update_lambda(
					[&](double i_lambda, double i_theta, double &io_data)
				{
						double d = std::acos(
								std::sin(theta_c)*std::sin(i_theta) +
								std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
						);

						double i_exp_fac = 20.0;
						io_data = std::exp(-d*d*i_exp_fac)*0.1*simVars->sim.h0;

						io_data *= simVars->sim.gravitation;
					}
				);

				/*
				 * Both versions are working
				 */
	#if 1
				SphereData stream_function(o_phi.sphereDataConfig);

				stream_function.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						double alpha = simVars->setup.advection_rotation_angle;

						io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
					}
				);

				o_vort = op->laplace(stream_function);

	#else

				o_vort.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						double alpha = simVars->setup.advection_rotation_angle;

						io_data = 2.0*u0/a*(-std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha) + std::sin(i_theta)*std::cos(alpha));
					}
				);
	#endif
				o_div.spectral_set_zero();
			}
			else if (
					simVars->setup.benchmark_name == "galewsky" ||			///< Standard Galewsky benchmark
					simVars->setup.benchmark_name == "galewsky_nobump" ||	///< Galewsky benchmark without bumps
					simVars->setup.benchmark_name == "galewsky_nosetparam"	///< Galewsky benchmark without overriding parameters
			)
			{

				if (simVars->setup.benchmark_name != "galewsky_nosetparam")
				{
					if (simVars->timecontrol.current_simulation_time == 0)
					{
						std::cout << "!!! WARNING !!!" << std::endl;
						std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
						std::cout << "!!! WARNING !!!" << std::endl;
					}

					/// Setup Galewski parameters
					simVars->sim.coriolis_omega = 7.292e-5;
					simVars->sim.gravitation = 9.80616;
					simVars->sim.earth_radius = 6.37122e6;
					simVars->sim.h0 = 10000;
				}

				op->setup(o_phi.sphereDataConfig, simVars->sim.earth_radius);

				/*
				 * Parameters from Galewsky paper setup
				 */
				double umax = 80.;

				double phi0 = M_PI/7.;
				double phi1 = 0.5*M_PI - phi0;
				double phi2 = 0.25*M_PI;		/// latitude placement of gaussian bump
				double en = std::exp(-4.0/std::pow((phi1-phi0), 2.0));
				double alpha = 1./3.;
				double beta = 1./15.;
				double hamp = 120.;

				if (simVars->setup.benchmark_galewsky_umax >= 0)
					umax = simVars->setup.benchmark_galewsky_umax;

				if (simVars->setup.benchmark_galewsky_hamp >= 0)
					hamp = simVars->setup.benchmark_galewsky_hamp;

				if (simVars->setup.benchmark_galewsky_phi2 >= 0)
					phi2 = simVars->setup.benchmark_galewsky_phi2;

				/*
				 * Setup V=0
				 */
				SphereDataPhysical vg(o_phi.sphereDataConfig);
				vg.physical_set_zero();

				/*
				 * Setup U=...
				 * initial velocity along longitude
				 */
				SphereDataPhysical ug(o_phi.sphereDataConfig);
				ug.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						if (phi >= phi1 || phi <= phi0)
							o_data = 0;
						else
							o_data = umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
					}
				);

				if (simVars->misc.sphere_use_robert_functions)
				{
					ug.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					vg.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					op->robert_uv_to_vortdiv(ug, vg, o_vort, o_div);
				}
				else
				{
					op->uv_to_vortdiv(ug, vg, o_vort, o_div);
				}

				computeGeostrophicBalance(
						o_phi,
						o_vort,
						o_div
				);

				SphereDataPhysical hbump(o_phi.sphereDataConfig);
				if (simVars->setup.benchmark_name == "galewsky")
				{
					hbump.physical_update_lambda(
						[&](double lon, double phi, double &o_data)
						{
							o_data = hamp*std::cos(phi)*std::exp(-(lon/alpha)*(lon/alpha))*std::exp(-(phi2-phi)*(phi2-phi)/beta);
						}
					);
				}

				o_phi += hbump*simVars->sim.gravitation;
			}
			else if (
					simVars->setup.benchmark_name == "williamson4"		||
					simVars->setup.benchmark_name == "forced_nonlinear"
			)
			{
				FatalError("TODO: Implement this");
			}
			else if (
					simVars->setup.benchmark_name == "williamson5"	||
					simVars->setup.benchmark_name == "flow_over_mountain"
			)
			{
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				/// Setup Williamson's parameters
				simVars->sim.coriolis_omega = 7.292e-5;
				simVars->sim.gravitation	= 9.80616;
				simVars->sim.earth_radius   = 6.37122e6;
				simVars->sim.h0			 = 5600;


				// update operator because we changed the simulation parameters
				op->setup(o_phi.sphereDataConfig, simVars->sim.earth_radius);

				const double u0 = 20.0;

				/*
				 * Setup V=0
				 */
				SphereDataPhysical vg(o_phi.sphereDataConfig);
				vg.physical_set_zero();

				/*
				 * Setup U=...
				 * initial velocity along longitude
				 */
				SphereDataPhysical ug(o_phi.sphereDataConfig);
				ug.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = u0 * std::cos(phi);
					}
				);

				if (simVars->misc.sphere_use_robert_functions)
				{
					ug.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					vg.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					op->robert_uv_to_vortdiv(ug, vg, o_vort, o_div);
				}
				else
				{
					op->uv_to_vortdiv(ug, vg, o_vort, o_div);
				}

				SphereDataPhysical hg(o_phi.sphereDataConfig);
				computeGeostrophicBalance(
						o_phi,
						o_vort,
						o_div
				);
			}
			else if (
					simVars->setup.benchmark_name == "williamson6"	||
					simVars->setup.benchmark_name == "rossby_haurwitz_wave"
			)
			{
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				/// Setup Williamson's parameters
				simVars->sim.coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.earth_radius = 6.37122e6;
				simVars->sim.h0 = 8000;

				// update operator because we changed the simulation parameters
				op->setup(o_phi.sphereDataConfig, simVars->sim.earth_radius);

				const double omega = 7.484e-6;
				const double K = omega;
				const int R	= 4;
				const double a = simVars->sim.earth_radius;

				/*
				 * Setup U=...
				 */
				SphereDataPhysical ug(o_phi.sphereDataConfig);
				ug.physical_update_lambda(
								[&](double lon, double phi, double &o_data)
									{
										o_data = a * omega * cos(phi)
												+ a * K * pow(cos(phi), R-1) * (R * sin(phi)*sin(phi) - cos(phi)*cos(phi)) * cos(R*lon);
									}
								);


				/*
				 * Setup V=...
				 */
				SphereDataPhysical vg(o_phi.sphereDataConfig);
				vg.physical_update_lambda(
								  [&](double lon, double phi, double &o_data)
										{
											o_data = - a * K * R * pow(cos(phi), R-1) * sin(phi) * sin(R*lon);
										}
								  );

				if (simVars->misc.sphere_use_robert_functions)
				{
					ug.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					vg.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					op->robert_uv_to_vortdiv(ug, vg, o_vort, o_div);
				}
				else
				{
					op->uv_to_vortdiv(ug, vg, o_vort, o_div);
				}

				SphereDataPhysical hg(o_phi.sphereDataConfig);
				computeGeostrophicBalance(
						o_phi,
						o_vort,
						o_div
				);
			}
			else if (
					simVars->setup.benchmark_name == "williamson7"	||
					simVars->setup.benchmark_name == "real_initial_conditions"
			)
			{
				FatalError("Williamson#7 not yet implemented!");
			}
			else if (simVars->setup.benchmark_name == "flat")
			{
				o_phi.physical_set_all_value(simVars->sim.h0*simVars->sim.gravitation);
				o_vort.physical_set_all_value(0);
				o_div.physical_set_all_value(0);
			}
			else if (simVars->setup.benchmark_name == "gaussian_bumps2")
			{
				SphereData tmp(o_phi.sphereDataConfig);
				SphereData o_h(o_phi.sphereDataConfig);

				o_h.physical_set_all_value(simVars->sim.h0);

				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, *simVars, 2.0*M_PI*0.1, M_PI/3, 20.0);
				o_h += (tmp-simVars->sim.h0);
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, *simVars, 2.0*M_PI*0.6, M_PI/5.0, 80.0);
				o_h += (tmp-simVars->sim.h0);
				BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, *simVars, 2.0*M_PI*0.8, -M_PI/4, 360.0);
				o_h += (tmp-simVars->sim.h0);

				o_phi = o_h*simVars->sim.gravitation;
			}
			else if (
					simVars->setup.benchmark_name == "geostrophic_balance"	||
					simVars->setup.benchmark_name == "geostrophic_balance_1"	||
					simVars->setup.benchmark_name == "geostrophic_balance_2"	||
					simVars->setup.benchmark_name == "geostrophic_balance_4"	||
					simVars->setup.benchmark_name == "geostrophic_balance_8"	||
					simVars->setup.benchmark_name == "geostrophic_balance_16"	||
					simVars->setup.benchmark_name == "geostrophic_balance_32"	||
					simVars->setup.benchmark_name == "geostrophic_balance_64"	||
					simVars->setup.benchmark_name == "geostrophic_balance_128"	||
					simVars->setup.benchmark_name == "geostrophic_balance_256"	||
					simVars->setup.benchmark_name == "geostrophic_balance_512"	||
					simVars->setup.benchmark_name == "geostrophic_balance_nosetparam"
			)
			{
				/*
				 * geostrophic_balance / geostrophic_balance_1:
				 * Williamson test case 2 for geostrophic balance.
				 *
				 * See Williamson paper for accurate setup
				 *
				 * "geostrophic_balance_N" means that N is the multiplier for the frequency
				 * in the direction of the Latitude
				 */
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				if (simVars->setup.benchmark_name != "geostrophic_balance_nosetparam")
				{
					simVars->sim.coriolis_omega = 7.292e-5;
					simVars->sim.gravitation = 9.80616;
					simVars->sim.earth_radius = 6.37122e6;
					simVars->sim.h0 = 29400.0/simVars->sim.gravitation;
				}

				SphereData o_h(o_phi.sphereDataConfig);

				// update operator because we changed the simulation parameters
				op->setup(o_h.sphereDataConfig, simVars->sim.earth_radius);

				//double u0 = (2.0*M_PI*simVars->sim.earth_radius)/(12.0*24.0*60.0*60.0);
				//double a = simVars->sim.earth_radius;

				//double phi0 = simVars->sim.h0*simVars->sim.gravitation;

				double freq_multiplier = 1.0;

				if (simVars->setup.benchmark_name == "geostrophic_balance_2")
					freq_multiplier = 2.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_4")
					freq_multiplier = 4.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_8")
					freq_multiplier = 8.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_16")
					freq_multiplier = 16.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_32")
					freq_multiplier = 32.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_64")
					freq_multiplier = 64.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_128")
					freq_multiplier = 128.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_256")
					freq_multiplier = 256.0;
				else if (simVars->setup.benchmark_name == "geostrophic_balance_512")
					freq_multiplier = 512.0;


				/*
				 * Setup V=0
				 */
				SphereDataPhysical vg(o_h.sphereDataConfig);
				vg.physical_set_zero();

				/*
				 * Setup U=...
				 * initial velocity along longitude
				 */
				SphereDataPhysical ug(o_h.sphereDataConfig);
				ug.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = std::cos(phi*freq_multiplier);
					}
				);

				if (simVars->misc.sphere_use_robert_functions)
				{
					ug.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					vg.physical_update_lambda_cosphi_grid(
						[&](double lon, double phi, double &o_data)
						{
							o_data *= phi;
						}
					);

					op->robert_uv_to_vortdiv(ug, vg, o_vort, o_div);
				}
				else
				{
					op->uv_to_vortdiv(ug, vg, o_vort, o_div);
				}

				computeGeostrophicBalance(
						o_phi,
						o_vort,
						o_div
				);
			}
			else
			{
				FatalError("Benchmark not implemented");
			}
		}
		else
		{
			SphereData h(o_phi.sphereDataConfig);
			SphereDataPhysical u(o_phi.sphereDataConfig);
			SphereDataPhysical v(o_phi.sphereDataConfig);

			setupInitialConditions_HUV(h, u, v);

			o_phi = h*simVars->sim.gravitation;
			if (simVars->misc.sphere_use_robert_functions)
				op->robert_uv_to_vortdiv(u, v, o_vort, o_div);
			else
				op->uv_to_vortdiv(u, v, o_vort, o_div);
		}
	}


public:
	void setupInitialConditions_HUV(
			SphereData &o_h,
			SphereDataPhysical &o_u,
			SphereDataPhysical &o_v
	)
	{
		if (simVars->setup.benchmark_id < 0)
		{
			printAvailableBenchmarks();
			FatalError("Benchmark scenario not selected");
		}

		if (simVars->setup.benchmark_id == 0)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, *simVars, M_PI/3.0, M_PI/3.0);
		}
		else if (simVars->setup.benchmark_id == 2)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, *simVars, M_PI/3.0, 0);
		}
		else if (simVars->setup.benchmark_id == 3)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, *simVars, M_PI/3.0, M_PI/3.0);
		}
		else if (simVars->setup.benchmark_id == 4)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, *simVars, M_PI/4.0, -M_PI);
		}
		else if (simVars->setup.benchmark_id == 5)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, *simVars, M_PI/4.0, -M_PI, 100.0);
		}
		else if (simVars->setup.benchmark_id == 6)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, *simVars, M_PI/3.0, M_PI/3.0, 20.0);
		}
		else if (simVars->setup.benchmark_id == 7)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, *simVars, M_PI/4.0, -M_PI, 100.0);
		}
		else if (simVars->setup.benchmark_id == 9)
		{
			SphereData tmp(o_h.sphereDataConfig);

			o_u.physical_set_all_value(0);
			o_v.physical_set_all_value(0);

			o_h.physical_set_all_value(simVars->sim.h0);

			BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, *simVars, 2.0*M_PI*0.1, M_PI/3, 20.0);
			o_h += (tmp-simVars->sim.h0);
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, *simVars, 2.0*M_PI*0.6, M_PI/5.0, 50.0);
			o_h += (tmp-simVars->sim.h0);
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, o_u, o_v, *simVars, 2.0*M_PI*0.8, -M_PI/4, 100.0);
			o_h += (tmp-simVars->sim.h0);
		}
		else if (simVars->setup.benchmark_id == 10 || simVars->setup.benchmark_id == 11)
		{
			/*
			 * Williamson test case 2 for geostrophic balance.
			 *
			 * See Williamson paper for accurate setup
			 */

			double u0 = 1.0;
			if (simVars->setup.benchmark_id == 11)
			{
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				simVars->sim.coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.earth_radius = 6.37122e6;
				simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

				op->setup(o_h.sphereDataConfig, simVars->sim.earth_radius);

				u0 = (2.0*M_PI*simVars->sim.earth_radius)/(12.0*24.0*60.0*60.0);
			}
			double a = simVars->sim.earth_radius;

			double phi0 = simVars->sim.h0*simVars->sim.gravitation;

			o_h.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
					io_data = (phi0 + (a*simVars->sim.coriolis_omega*u0*std::cos(i_lat)*std::cos(i_lat)))/simVars->sim.gravitation;
				}
			);

			if (simVars->misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						io_data = u0*std::cos(i_lat)*std::cos(i_lat);
					}
				);
			}
			else
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						io_data = std::cos(i_lat);
					}
				);
			}

			o_v.physical_set_zero();

		}
		else if (simVars->setup.benchmark_id == 21)
		{
			/*
			 * Advection test case
			 * See Williamson test case, eq. (75), (76)
			 *
			 * phi_t = DIV (u phi)
			 */

			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			simVars->sim.coriolis_omega = 7.292e-5;
			simVars->sim.gravitation = 9.80616;
			simVars->sim.earth_radius = 6.37122e6;
			simVars->sim.h0 = 1000.0;

			op->setup(o_h.sphereDataConfig, simVars->sim.earth_radius);

			double lambda_c = 3.0*M_PI/2.0;
			double theta_c = M_PI/4.0;
			double a = simVars->sim.earth_radius;

			double R = a/3.0;
			double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

			o_h.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
			{
					double r = a * std::acos(
							std::sin(theta_c)*std::sin(i_theta) +
							std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
					);

					if (r < R)
						io_data = simVars->sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
					else
						io_data = 0;
				}
			);


			if (simVars->misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
							u0*(
									std::cos(i_theta)*std::cos(simVars->setup.advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
								);

						io_data *= std::cos(i_lat);
					}
				);

				o_v.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
//						double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
									std::sin(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
								);

						io_data *= std::cos(i_lat);
					}
				);
			}
			else
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
							u0*(
								std::cos(i_theta)*std::cos(simVars->setup.advection_rotation_angle) +
								std::sin(i_theta)*std::cos(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
								);
					}
				);

				o_v.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
//						double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
								std::sin(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
							);
					}
				);
			}

			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Storing advection in output velocities !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			std::cout << "advection_rotation_angle: " << simVars->setup.advection_rotation_angle << std::endl;
		}
		else if (simVars->setup.benchmark_id == 22)
		{
			/**
			 * Advection test case
			 * See Williamson test case, eq. (75), (76)
			 *
			 * phi_t = phi GRAD u
			 */

			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			simVars->sim.coriolis_omega = 7.292e-5;
			simVars->sim.gravitation = 9.80616;
			simVars->sim.earth_radius = 6.37122e6;
			simVars->sim.h0 = 1000.0;

			op->setup(o_h.sphereDataConfig, simVars->sim.earth_radius);

			double lambda_c = 3.0*M_PI/2.0;
			double theta_c = 0;
			double a = simVars->sim.earth_radius;

			double R = a/3.0;
			double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

			o_h.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					double r = a * std::acos(
							std::sin(theta_c)*std::sin(i_theta) +
							std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
					);

					if (r < R)
						io_data = simVars->sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
					else
						io_data = 0;
				}
			);


			if (simVars->misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
								u0*(
									std::cos(i_theta)*std::cos(simVars->setup.advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
							);

						io_data *= std::cos(i_lat);
					}
				);

				o_v.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
//						double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
									std::sin(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
							);

						io_data *= std::cos(i_lat);
					}
				);
			}
			else
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
								u0*(
									std::cos(i_theta)*std::cos(simVars->setup.advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
							);
					}
				);

				o_v.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
//						double i_phi = i_lat;
						double i_lambda = i_lon;
						io_data =
							-u0*(
									std::sin(i_lambda)*std::sin(simVars->setup.advection_rotation_angle)
							);
					}
				);
			}

			//o_h.physical_truncate();
			//o_u.physical_truncate();
			//o_v.physical_truncate();

			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Storing advection in output velocities !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}


			std::cout << "advection_rotation_angle: " << simVars->setup.advection_rotation_angle << std::endl;
		}
		else if (simVars->setup.benchmark_id == 40)
		{
			o_h.physical_set_all_value(simVars->sim.h0);
			o_u.physical_set_all_value(0);
			o_v.physical_set_all_value(0);
		}
		else if (simVars->setup.benchmark_id == 50)
		{
			/*
			 * PAUL N. SWARZTRAUBER
			 * Shallow Water Flow on the Sphere
			 * 2004
			 */
			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			if (simVars->sim.coriolis_omega != 0)
				simVars->sim.coriolis_omega = 7.292e-5;

			simVars->sim.gravitation = 9.80616;
			simVars->sim.earth_radius = 6.37122e6;
			simVars->sim.h0 = 29400.0;

			op->setup(o_h.sphereDataConfig, simVars->sim.earth_radius);

			o_u.physical_set_zero();
			o_v.physical_set_zero();

			const double a = simVars->sim.earth_radius;
			const double A = 6000.0;
			const double alpha = 10;

			const double center_lat = M_PI/4;
			const double center_lon = M_PI;

			o_h.physical_update_lambda(
				[&](double i_lambda, double i_phi, double &io_data)
				{
				  double x = a* ( std::cos(i_phi)*std::cos(i_lambda) - std::cos(center_lat)*std::cos(center_lon) );
				  double y = a* ( std::cos(i_phi)*std::sin(i_lambda) - std::cos(center_lat)*std::sin(center_lon) );
				  double z = a* ( std::sin(i_phi)					- std::sin(center_lat) );

					double d = std::sqrt(x*x+y*y+z*z);

					io_data = simVars->sim.h0 + A*std::exp(-alpha*(d/a)*(d/a));
				}
			);
		}
		else if (simVars->setup.benchmark_id == 200)
		{
			o_h.physical_set_all_value(simVars->sim.h0);
			o_u.physical_set_zero();
			o_v.physical_set_zero();
		}
		else
		{
			FatalError("Unknown scenario id");
		}
	}
};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SWESPHEREBENCHMARKSCOMBINED_HPP_ */
