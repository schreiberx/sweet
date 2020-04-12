/*
 * SWESphereBenchmarksCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *	  Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SPHERE_SWESPHEREBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_SPHERE_SWESPHEREBENCHMARKSCOMBINED_HPP_

#include <functional>
#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/BenchmarkGaussianDam.hpp>
#include <benchmarks_sphere/BenchmarkFlowOverMountain.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <libmath/GaussQuadrature.hpp>

#if SWEET_MPI
	#include <mpi.h>
#endif

#if !SWEET_USE_SPHERE_SPECTRAL_SPACE
	#error "SWEET_USE_SPHERE_SPECTRAL_SPACE not enabled"
#endif

class SWESphereBenchmarksCombined
{
public:
	// Simulation variables
	SimulationVariables *simVars;

	// Sphere operators
	SphereOperators_SphereData *op;

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
		std::cout << "         OPTION:" << std::endl;
		std::cout << "         --advection-rotation-angle=[angle]" << std::endl;
		std::cout << std::endl;
		std::cout << "  WILLIAMSON #1 (alternative):" << std::endl;
		std::cout << "     'adv_gauss_bump': Advection test case of gaussian bump" << std::endl;
		std::cout << "         OPTION:" << std::endl;
		std::cout << "         --advection-rotation-angle=[angle]" << std::endl;
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



	void normalization(
			SphereData_Spectral &io_phi,
			const std::string &normalization = ""
	)
	{
		if (normalization == "avg_zero")
		{
			// move average value to 0
			double phi_min = io_phi.getSphereDataPhysical().physical_reduce_min();
			double phi_max = io_phi.getSphereDataPhysical().physical_reduce_max();

			double avg = 0.5*(phi_max+phi_min);

			io_phi -= avg;
		}
		else if (normalization == "min_zero")
		{
			// move minimum value to zero
			double phi_min = io_phi.getSphereDataPhysical().physical_reduce_min();
			io_phi -= phi_min;
		}
		else if (normalization == "max_zero")
		{
			// move maximum value to zero
			double phi_max = io_phi.getSphereDataPhysical().physical_reduce_max();
			io_phi -= phi_max;
		}
		else if (normalization == "")
		{
		}
		else
		{
			FatalError("Normalization not supported!");
		}
	}



	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_nonlinear(
			SphereData_Spectral &i_vort,
			SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi
	)
	{
		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereData_Physical ug(o_phi.sphereDataConfig);
		SphereData_Physical vg(o_phi.sphereDataConfig);

		op->vortdiv_to_uv(i_vort, i_div, ug, vg);

		SphereData_Physical vrtg = i_vort.getSphereDataPhysical();

		SphereData_Physical tmpg1 = ug*(vrtg+op->fg);
		SphereData_Physical tmpg2 = vg*(vrtg+op->fg);

		SphereData_Spectral tmpspec1(o_phi.sphereDataConfig);
		SphereData_Spectral tmpspec2(o_phi.sphereDataConfig);

		op->uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		SphereData_Spectral phispec = op->inv_laplace(tmpspec1) - 0.5*(ug*ug+vg*vg);

		o_phi = phispec;
	}



	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_linear(
			SphereData_Spectral &i_vort,
			SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi
	)
	{
		const SphereData_Config *sphereDataConfig = o_phi.sphereDataConfig;
		/*
		 * Setup Coriolis effect
		 */
		SphereData_Physical f(sphereDataConfig);
		f.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = 2.0*simVars->sim.sphere_rotating_coriolis_omega*mu;
			}
		);

		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereData_Physical u(sphereDataConfig);
		SphereData_Physical v(sphereDataConfig);

		op->vortdiv_to_uv(i_vort, i_div, u, v);

		SphereData_Physical vrtg = i_vort.getSphereDataPhysical();

		SphereData_Physical tmpg1 = u*f;
		SphereData_Physical tmpg2 = v*f;

		SphereData_Spectral tmpspec1(sphereDataConfig);
		SphereData_Spectral tmpspec2(sphereDataConfig);

		op->uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		SphereData_Spectral phispec = op->inv_laplace(tmpspec1);

		o_phi = simVars->sim.h0*simVars->sim.gravitation + phispec;
	}




public:
	void setupTopography()
	{
		if (simVars == nullptr)
			FatalError("Benchmarks are not yet initialized");

		if (simVars->benchmark.benchmark_name == "flow_over_mountain")
		{
			// set the topography flag to true
			simVars->benchmark.use_topography = true;

			// setup the parameters for the flow-over-mountain test case
			const double R			= M_PI/9.;
			const double h_topo_0	 = 2000.;
			const double i_center_lon = 3.*M_PI/2.;
			const double i_center_lat = M_PI/6.;

			// initialize the topography
			simVars->benchmark.h_topo.physical_set_zero();

			// setup the topography vector
			BenchmarkFlowOverMountain::setup_topography(
					simVars->benchmark.h_topo,
					*simVars,
					R,
					h_topo_0,
					i_center_lon,
					i_center_lat
			);
		}
		else
		{
			// set the topography flag to false
			simVars->benchmark.use_topography = false;
		}
	}



	void setup(
			SimulationVariables &io_simVars,
			SphereOperators_SphereData &io_op
	)
	{
		simVars = &io_simVars;
		op = &io_op;
	}

	void setupInitialConditions_pert(
			SphereData_Spectral &o_phi_pert,
			SphereData_Spectral &o_vort,
			SphereData_Spectral &o_div
	)
	{
		double gh0 = simVars->sim.gravitation*simVars->sim.h0;

		if (simVars == nullptr)
			FatalError("Benchmarks are not yet initialized");

		simVars->iodata.output_time_scale = 1.0/(60.0*60.0);

		if (simVars->benchmark.benchmark_name == "")
			FatalError("Benchmark name not specified!");

		const SphereData_Config *sphereDataConfig = o_phi_pert.sphereDataConfig;

		if (simVars->benchmark.benchmark_name == "gaussian_bump_advection")
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
				SphereData_Spectral* o_sphere_data = (SphereData_Spectral*)o_data_void;
				SWESphereBenchmarksCombined* s = (SWESphereBenchmarksCombined*)o_data_user_void;

				if (i_field_id == 1)
				{
					double a = s->simVars->sim.sphere_radius;
					double alpha = s->simVars->benchmark.sphere_advection_rotation_angle;
					double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

					double r;
					if (s->simVars->sim.advection_velocity[2] == 0)
						r = 0;
					else
						r = i_simulation_timestamp/s->simVars->sim.advection_velocity[2]*2.0*M_PI;

					// change velocity
					u0 = u0*(1.0 + std::cos(r));

					SphereData_Physical stream_function(o_sphere_data->sphereDataConfig);

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

			simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
			simVars->sim.gravitation = 9.80616;
			simVars->sim.sphere_radius = 6.37122e6;
			simVars->sim.h0 = 1000.0;

			op->setup(sphereDataConfig, &(simVars->sim));

			double lambda_c = 3.0*M_PI/2.0;
			double theta_c = 0.0;

			SphereData_Physical phi_phys(o_phi_pert.sphereDataConfig);

			phi_phys.physical_update_lambda(
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

			o_phi_pert.loadSphereDataPhysical(phi_phys - gh0);


			if (simVars->sim.advection_velocity[2] != 0)
			{
				// backup data config
				ext_forces_data_config = o_phi_pert.sphereDataConfig;

				// set callback
				simVars->benchmark.getExternalForcesCallback = callback_external_forces_advection_field;

				// set user data to this class
				simVars->benchmark.getExternalForcesUserData = this;
			}

			// setup velocities with initial time stamp
			callback_external_forces_advection_field(1, simVars->timecontrol.current_simulation_time, &o_vort, this);
			callback_external_forces_advection_field(2, simVars->timecontrol.current_simulation_time, &o_div, this);
		}
		else if (
				simVars->benchmark.benchmark_name == "williamson1"		||
				simVars->benchmark.benchmark_name == "adv_cosine_bell"
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

			simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
			simVars->sim.gravitation = 9.80616;
			simVars->sim.sphere_radius = 6.37122e6;
			simVars->sim.h0 = 1000.0;

			// reset operator
			op->setup(o_phi_pert.sphereDataConfig, &(simVars->sim));

			double lambda_c = 3.0*M_PI/2.0;
			double theta_c = 0.0;
			double a = simVars->sim.sphere_radius;

			double R = a/3.0;
			double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

			SphereData_Physical phi_phys(o_phi_pert.sphereDataConfig);

			phi_phys.physical_update_lambda(
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

			o_phi_pert.loadSphereDataPhysical(phi_phys);

			SphereData_Physical stream_function(o_phi_pert.sphereDataConfig);

			stream_function.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
					double i_theta = i_lat;
					double i_lambda = i_lon;
					double alpha = simVars->benchmark.sphere_advection_rotation_angle;

					io_data = -a*u0*(std::sin(i_theta)*std::cos(alpha) - std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha));
				}
			);

			o_vort = op->laplace(stream_function);
			o_div.spectral_set_zero();

			o_phi_pert -= gh0;
		}
		else if (
				simVars->benchmark.benchmark_name == "williamson1b"		||
				simVars->benchmark.benchmark_name == "adv_gauss_bump"
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

			simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
			simVars->sim.gravitation = 9.80616;
			simVars->sim.sphere_radius = 6.37122e6;
			simVars->sim.h0 = 1000.0;

			op->setup(o_phi_pert.sphereDataConfig, &(simVars->sim));

			double lambda_c = 3.0*M_PI/2.0;
			double theta_c = 0.0;
			double a = simVars->sim.sphere_radius;

			//double R = a/3.0;
			double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);
			//double u0 = (2.0*M_PI*a*1000.0);

			SphereData_Physical phi_phys(o_phi_pert.sphereDataConfig);

			phi_phys.physical_update_lambda(
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

			o_phi_pert.loadSphereDataPhysical(phi_phys);

			/*
			 * Both versions are working
			 */
#if 1
			SphereData_Physical stream_function(o_phi_pert.sphereDataConfig);

			stream_function.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
					double i_theta = i_lat;
					double i_lambda = i_lon;
					double alpha = simVars->benchmark.sphere_advection_rotation_angle;

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
					double alpha = simVars->benchmark.sphere_advection_rotation_angle;

					io_data = 2.0*u0/a*(-std::cos(i_lambda)*std::cos(i_theta)*std::sin(alpha) + std::sin(i_theta)*std::cos(alpha));
				}
			);
#endif
			o_div.spectral_set_zero();

			o_phi_pert -= gh0;
		}
		else if (
				simVars->benchmark.benchmark_name == "galewsky" ||			///< Standard Galewsky benchmark
				simVars->benchmark.benchmark_name == "galewsky_nobump" ||	///< Galewsky benchmark without bumps
				simVars->benchmark.benchmark_name == "galewsky_nosetparams"	///< Galewsky benchmark without overriding parameters
		)
		{
			if (simVars->benchmark.benchmark_name != "galewsky_nosetparams")
			{
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				/// Setup Galewski parameters
				simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.sphere_radius = 6.37122e6;

				// see doc/galewsky_mean_layer_depth/ on how to get this constant.
				// it is NOT 10e3 (see Galewsky paper)
				simVars->sim.h0 = 10158.186170454619;


				/*
				 * Rerun setup to update the operators with the potentially new values
				 */
				op->setup(o_phi_pert.sphereDataConfig, &(simVars->sim));
			}

			/*
			 * Parameters from Galewsky paper setup
			 */
			double a = simVars->sim.sphere_radius;
			double omega = simVars->sim.sphere_rotating_coriolis_omega;
			double umax = 80.;
			double phi0 = M_PI/7.;
			double phi1 = 0.5*M_PI - phi0;
			double phi2 = 0.25*M_PI;		/// latitude placement of gaussian bump
			double en = std::exp(-4.0/std::pow((phi1-phi0), 2.0));
			double alpha = 1./3.;
			double beta = 1./15.;
			double hamp = 120.;


			if (simVars->benchmark.benchmark_galewsky_umax >= 0)
				umax = simVars->benchmark.benchmark_galewsky_umax;

			if (simVars->benchmark.benchmark_galewsky_hamp >= 0)
				hamp = simVars->benchmark.benchmark_galewsky_hamp;

			if (simVars->benchmark.benchmark_galewsky_phi2 >= 0)
				phi2 = simVars->benchmark.benchmark_galewsky_phi2;

			/*
			 * Setup V=0
			 */
			SphereData_Physical vg(o_phi_pert.sphereDataConfig);
			vg.physical_set_zero();

			auto lambda_u = [&](double phi) -> double
			{
				if (phi >= phi1-1e-5 || phi <= phi0+1e-5)
					return 0.0;
				else
					return umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
			};

			auto lambda_f = [&](double phi) -> double
			{
				return a*lambda_u(phi)*(2.0*omega*std::sin(phi)+(std::tan(phi)/a)*lambda_u(phi));
			};

			/*
			 * Setup U=...
			 * initial velocity along longitude
			 */
			SphereData_Physical ug(o_phi_pert.sphereDataConfig);
			ug.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					o_data = lambda_u(phi);
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

			bool use_analytical_geostrophic_setup = simVars->misc.comma_separated_tags.find("galewsky_analytical_geostrophic_setup") != std::string::npos;

			if (use_analytical_geostrophic_setup)
			{
				std::cout << "[MULE] use_analytical_geostrophic_setup: 1" << std::endl;
				computeGeostrophicBalance_nonlinear(
						o_vort,
						o_div,
						o_phi_pert
				);

				double h0_ = 10e3;
				o_phi_pert = simVars->sim.gravitation * h0_ + o_phi_pert;


			}
			else
			{
				std::cout << "[MULE] use_analytical_geostrophic_setup: 0" << std::endl;

				/*
				 * Initialization of SWE height
				 *
				 * Metric correction terms based on John Thuburn's code
				 */
#if 1
				const unsigned short nlat = sphereDataConfig->physical_num_lat;
				std::vector<double> hg_cached;
				hg_cached.resize(nlat);

				double h_metric_area = 0;
				double hg_sum = 0;
				double int_start, int_end, int_delta;

				int j = sphereDataConfig->physical_num_lat-1;


				// start/end of first integration interval
				{
					assert(sphereDataConfig->lat[j] < 0);

					// start at the south pole
					int_start = -M_PI*0.5;

					// first latitude gaussian point
					int_end = sphereDataConfig->lat[j];

					// 1d area of integration
					int_delta = int_end - int_start;

					assert(int_delta > 0);
					assert(int_delta < 1);

					double hg = GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 20);
					//hg = (int_end+int_start)*0.5;
					hg_cached[j] = hg;

					/*
					 * cos scaling is required for 2D sphere coverage at this latitude
					 *
					 * metric term which computes the area coverage of each point
					 */
					// use integrated average as below instead of the following formulation
					// double mterm = cos((int_start+int_end)*0.5);
					//double mterm = (std::sin(int_end)-std::sin(int_start))*2.0*M_PI;
					double mterm = std::cos(sphereDataConfig->lat[j])*2.0*M_PI;
					assert(mterm > 0);

					hg_sum += hg*mterm;
					h_metric_area += mterm;

					int_start = int_end;
				}
				j--;

				for (; j >= 0; j--)
				{
					double int_end = sphereDataConfig->lat[j];
					int_delta = int_end - int_start;
					assert(int_delta > 0);

					double hg = hg_cached[j+1] + GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 20);

					//hg = (int_end+int_start)*0.5;
					hg_cached[j] = hg;

					// metric term which computes the area coverage of each point
					//double mterm = (std::sin(int_end)-std::sin(int_start))*2.0*M_PI;
					double mterm = std::cos(sphereDataConfig->lat[j])*2.0*M_PI;

					hg_sum += hg*mterm;
					h_metric_area += mterm;

					// continue at the end of the last integration interval
					int_start = int_end;
				}

				// last integration interval
				{
					assert(int_start > 0);
					int_end = M_PI*0.5;

					int_delta = int_end - int_start;
					assert(int_delta > 0);

					// metric term which computes the area coverage of each point
					//double mterm = (std::sin(int_end)-std::sin(int_start))*2.0*M_PI;
					double mterm = std::cos(sphereDataConfig->lat[0])*2.0*M_PI;

					double hg = hg_cached[0] + GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 20);
					//hg = (int_end+int_start)*0.5;
					hg_sum += hg*mterm;
					h_metric_area += mterm;
				}

				assert(h_metric_area > 0);

				double h_sum = hg_sum / simVars->sim.gravitation;
				double h_comp_avg = h_sum / h_metric_area;

				double h0 = 10000.0 + h_comp_avg;
				std::cout << "Galewsky benchmark H0 (computed, not used!): " << h0 << std::endl;

#else

				std::vector<double> hg_cached;
				hg_cached.resize(sphereDataConfig->physical_num_lat);

				double int_start = -M_PI*0.5;
				for (int j = sphereDataConfig->physical_num_lat-1; j >= 0; j--)
				{
					double int_end = sphereDataConfig->lat[j];
					double quad = GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 5);

					if (j == sphereDataConfig->physical_num_lat-1)
						hg_cached[j] = quad;
					else
						hg_cached[j] = hg_cached[j+1] + quad;

					int_start = int_end;


					std::cout << sphereDataConfig->lat[j] << ": " << hg_cached[j] << std::endl;
				}

#endif

				// update data
				SphereData_Physical phig(sphereDataConfig);
				phig.physical_update_lambda_array(
					[&](int i, int j, double &o_data)
					{
						o_data = hg_cached[j];
					}
				);

				o_phi_pert.loadSphereDataPhysical(phig);

				o_phi_pert = gh0 - o_phi_pert;
			}

#if 1


			SphereData_Physical hbump(o_phi_pert.sphereDataConfig);
			if (simVars->benchmark.benchmark_name == "galewsky")
			{
				hbump.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = hamp*std::cos(phi)*std::exp(-std::pow((lon-M_PI)/alpha, 2.0))*std::exp(-std::pow((phi2-phi)/beta, 2.0));
					}
				);
				o_phi_pert += hbump*simVars->sim.gravitation;
			}
#endif

			o_phi_pert -= gh0;


			std::cout << gh0 << std::endl;
			std::cout << "phi min: " << o_phi_pert.getSphereDataPhysical().physical_reduce_min() << std::endl;
			std::cout << "phi max: " << o_phi_pert.getSphereDataPhysical().physical_reduce_max() << std::endl;
		}
		else if (
				simVars->benchmark.benchmark_name == "williamson4"		||
				simVars->benchmark.benchmark_name == "forced_nonlinear"
		)
		{
			FatalError("TODO: Implement this");
		}
		else if (
				simVars->benchmark.benchmark_name == "williamson5"	||
				simVars->benchmark.benchmark_name == "flow_over_mountain"
		)
		{
			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			/// Setup Williamson's parameters
			simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
			simVars->sim.gravitation	= 9.80616;
			simVars->sim.sphere_radius   = 6.37122e6;
			simVars->sim.h0			 = 5600;


			// update operator because we changed the simulation parameters
			op->setup(o_phi_pert.sphereDataConfig, &(simVars->sim));

			const double u0 = 20.0;

			/*
			 * Setup V=0
			 */
			SphereData_Physical vg(o_phi_pert.sphereDataConfig);
			vg.physical_set_zero();

			/*
			 * Setup U=...
			 * initial velocity along longitude
			 */
			SphereData_Physical ug(o_phi_pert.sphereDataConfig);
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

			SphereData_Physical hg(o_phi_pert.sphereDataConfig);
			computeGeostrophicBalance_nonlinear(
					o_vort,
					o_div,
					o_phi_pert
			);

			o_phi_pert -= gh0;
		}
		else if (
				simVars->benchmark.benchmark_name == "williamson6"	||
				simVars->benchmark.benchmark_name == "rossby_haurwitz_wave"
		)
		{
			if (simVars->timecontrol.current_simulation_time == 0)
			{
				std::cout << "!!! WARNING !!!" << std::endl;
				std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
				std::cout << "!!! WARNING !!!" << std::endl;
			}

			/// Setup Williamson's parameters
			simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
			simVars->sim.gravitation = 9.80616;
			simVars->sim.sphere_radius = 6.37122e6;
			simVars->sim.h0 = 8000;

			// update operator because we changed the simulation parameters
			op->setup(o_phi_pert.sphereDataConfig, &(simVars->sim));

			const double omega = 7.484e-6;
			const double K = omega;
			const int R	= 4;
			const double a = simVars->sim.sphere_radius;

			/*
			 * Setup U=...
			 */
			SphereData_Physical ug(o_phi_pert.sphereDataConfig);
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
			SphereData_Physical vg(o_phi_pert.sphereDataConfig);
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

			SphereData_Physical hg(o_phi_pert.sphereDataConfig);
			computeGeostrophicBalance_nonlinear(
					o_vort,
					o_div,
					o_phi_pert
			);

			o_phi_pert -= gh0;
		}
		else if (
				simVars->benchmark.benchmark_name == "williamson7"	||
				simVars->benchmark.benchmark_name == "real_initial_conditions"
		)
		{
			FatalError("Williamson#7 not yet implemented!");
		}
		else if (simVars->benchmark.benchmark_name == "flat")
		{
			o_phi_pert.spectral_set_zero();
			o_phi_pert += simVars->sim.h0*simVars->sim.gravitation;
			o_vort.spectral_set_zero();
			o_div.spectral_set_zero();

			o_phi_pert -= gh0;
		}
		else if (simVars->benchmark.benchmark_name == "gaussian_bumps2" || simVars->benchmark.benchmark_name == "three_gaussian_bumps")
		{
			SphereData_Physical tmp(o_phi_pert.sphereDataConfig);
			SphereData_Spectral o_h(o_phi_pert.sphereDataConfig);

			tmp.physical_set_all_value(simVars->sim.h0);
			o_h.loadSphereDataPhysical(tmp);

			BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, *simVars, 2.0*M_PI*0.1, M_PI/3, 20.0);
			o_h += (SphereData_Spectral(tmp)-simVars->sim.h0);
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, *simVars, 2.0*M_PI*0.6, M_PI/5.0, 80.0);
			o_h += (SphereData_Spectral(tmp)-simVars->sim.h0);
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(tmp, *simVars, 2.0*M_PI*0.8, -M_PI/4, 360.0);
			o_h += (SphereData_Spectral(tmp)-simVars->sim.h0);

			o_phi_pert = o_h*simVars->sim.gravitation;

			o_phi_pert -= gh0;

		}
		else if (simVars->benchmark.benchmark_name == "gaussian_bumps_phi_vort_div")
		{
			{
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.sphere_radius = 6.37122e6;
				simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

				// Scale geopotential to make NL influencing the stiffness stronger
				simVars->sim.h0 *= 0.2;
				simVars->sim.gravitation *= 0.2;

				op->setup(sphereDataConfig, &(simVars->sim));
			}

			SphereData_Physical tmp(o_phi_pert.sphereDataConfig);

			BenchmarkGaussianDam::setup_initial_conditions_gaussian_normalized(tmp, *simVars, 2.0*M_PI*0.1, M_PI/3, 1.0);
			o_phi_pert.loadSphereDataPhysical(tmp);
			o_phi_pert *= 0.1;
			o_phi_pert += simVars->sim.h0*simVars->sim.gravitation;

			BenchmarkGaussianDam::setup_initial_conditions_gaussian_normalized(tmp, *simVars, 2.0*M_PI*0.1, M_PI/3, 1.0);
			o_vort.loadSphereDataPhysical(tmp);
			o_vort *= -1e-8;
			//o_vort *= 0;
			BenchmarkGaussianDam::setup_initial_conditions_gaussian_normalized(tmp, *simVars, 2.0*M_PI*0.1, M_PI/3, 1.0);
			o_div.loadSphereDataPhysical(tmp);
			o_div *= 1e-8;

			/*
			 * Convert forward/backward to velocity space to apply a certain truncation
			 */
			SphereData_Physical ug(o_phi_pert.sphereDataConfig);
			SphereData_Physical vg(o_phi_pert.sphereDataConfig);
			if (simVars->misc.sphere_use_robert_functions)
				op->robert_vortdiv_to_uv(o_vort, o_div, ug, vg);
			else
				op->vortdiv_to_uv(o_vort, o_div, ug, vg);
			if (simVars->misc.sphere_use_robert_functions)
				op->robert_uv_to_vortdiv(ug, vg, o_vort, o_div);
			else
				op->uv_to_vortdiv(ug, vg, o_vort, o_div);

			o_phi_pert -= gh0;
		}
		else if (
				simVars->benchmark.benchmark_name == "williamson2"			||
				simVars->benchmark.benchmark_name == "geostrophic_balance"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_1"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_2"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_4"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_8"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_16"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_32"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_64"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_128"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_256"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_512"	||
				simVars->benchmark.benchmark_name == "geostrophic_balance_nosetparams"
		)
		{
			/*
			 * geostrophic_balance / geostrophic_balance_1:
			 * Williamson test case 2 for geostrophic balance.
			 *
			 * WARNING: This test uses a balanced solution for the full non-linear equations
			 * See Williamson paper for accurate setup
			 *
			 * "geostrophic_balance_N" means that N is the multiplier for the frequency
			 * in the direction of the Latitude
			 */
			if (simVars->benchmark.benchmark_name != "geostrophic_balance_nosetparams")
			{
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.sphere_radius = 6.37122e6;
				simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

				op->setup(sphereDataConfig, &(simVars->sim));
			}

			double a = simVars->sim.sphere_radius;
			double omega = simVars->sim.sphere_rotating_coriolis_omega;
			double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
			double alpha = 0;


			double freq_multiplier = 1.0;

			if (simVars->benchmark.benchmark_name == "geostrophic_balance_2")
				freq_multiplier = 2.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_4")
				freq_multiplier = 4.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_8")
				freq_multiplier = 8.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_16")
				freq_multiplier = 16.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_32")
				freq_multiplier = 32.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_64")
				freq_multiplier = 64.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_128")
				freq_multiplier = 128.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_256")
				freq_multiplier = 256.0;
			else if (simVars->benchmark.benchmark_name == "geostrophic_balance_512")
				freq_multiplier = 512.0;


			/*
			 * Setup U
			 */
			SphereData_Physical ug(sphereDataConfig);
			ug.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 90, Williamson TC paper
					o_data = u0*(std::cos(phi*freq_multiplier)*std::cos(alpha) + std::cos(lon)*std::sin(phi*freq_multiplier)*std::sin(alpha));
				}
			);

			/*
			 * Setup V
			 */
			SphereData_Physical vg(sphereDataConfig);
			vg.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 91, Williamson TC paper
					o_data = -u0*std::sin(lon*freq_multiplier)*std::sin(alpha);
				}
			);

			op->uv_to_vortdiv(ug, vg, o_vort, o_div);



			/**
			 * TEST for non-divergent test case
			 */
			double div_zero = o_div.getSphereDataPhysical().physical_reduce_max_abs();
			if (div_zero > 1e-12)
			{

				std::cout << "Divergence: " << div_zero << std::endl;
				FatalError("Divergence should be close to 0, maybe there are some numerical round-off errors?");
			}

			/**
			 * TEST for correct vorticity
			 */
			/*
			 * Setup relative vorticity
			 */
			SphereData_Physical vortg(sphereDataConfig);
			vortg.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 94, Williamson TC paper

					// absolute vorticity, but we like the relative one
					//o_data = (2.0*u0/a + 2.0*omega)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));

					// relative vorticity
					o_data = (2.0*u0/a)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));
				}
			);

			double vort_diff = (o_vort.getSphereDataPhysical() - vortg).physical_reduce_max_abs();
			if (vort_diff > 1e-12)
			{

				std::cout << "Vorticity difference: " << vort_diff << std::endl;
				FatalError("Vorticity fields differ (should be close to 0), maybe there are some numerical round-off errors?");
			}


			bool use_analytical_geostrophic_setup = simVars->misc.comma_separated_tags.find("geostrophic_balance_analytical_setup") != std::string::npos;


			if (use_analytical_geostrophic_setup)
			{
				std::cout << "[MULE] geostrophic_balance_analytical_setup: 1" << std::endl;

				computeGeostrophicBalance_nonlinear(
						o_vort,
						o_div,
						o_phi_pert
				);
			}
			else
			{
				std::cout << "[MULE] geostrophic_balance_analytical_setup: 0" << std::endl;

				// Squared term in Eq. 95, Williamson TC paper
				SphereData_Physical r2(sphereDataConfig);
				r2.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = -std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha);
						o_data = o_data*o_data;
					}
				);

				// Eq. 95, Williamson TC paper
				SphereData_Physical phig = -(a*omega*u0 + u0*u0/2.0)*r2;

				o_phi_pert.loadSphereDataPhysical(phig);
			}
		}

		else if (
				simVars->benchmark.benchmark_name == "williamson2_linear"			||
				simVars->benchmark.benchmark_name == "geostrophic_balance_linear"
		)
		{
			/*
			 * Linear-only!!!
			 *
			 * The original version is for the non-linear case
			 */
			if (simVars->benchmark.benchmark_name != "geostrophic_balance_nosetparam")
			{
				if (simVars->timecontrol.current_simulation_time == 0)
				{
					std::cout << "!!! WARNING !!!" << std::endl;
					std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
					std::cout << "!!! WARNING !!!" << std::endl;
				}

				simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
				simVars->sim.gravitation = 9.80616;
				simVars->sim.sphere_radius = 6.37122e6;
				simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

				op->setup(sphereDataConfig, &(simVars->sim));
			}

			double a = simVars->sim.sphere_radius;
			//double omega = simVars->sim.sphere_rotating_coriolis_omega;
			double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
			double alpha = 0;

			double freq_multiplier = 1.0;


			/*
			 * Setup U
			 */
			SphereData_Physical ug(sphereDataConfig);
			ug.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 90, Williamson TC paper
					o_data = u0*(std::cos(phi*freq_multiplier)*std::cos(alpha) + std::cos(lon)*std::sin(phi*freq_multiplier)*std::sin(alpha));
				}
			);

			/*
			 * Setup V
			 */
			SphereData_Physical vg(sphereDataConfig);
			vg.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 91, Williamson TC paper
					o_data = -u0*std::sin(lon*freq_multiplier)*std::sin(alpha);
				}
			);

			op->uv_to_vortdiv(ug, vg, o_vort, o_div);



			/**
			 * TEST for non-divergent test case
			 */
			double div_zero = o_div.getSphereDataPhysical().physical_reduce_max_abs();
			if (div_zero > 1e-12)
			{

				std::cout << "Divergence: " << div_zero << std::endl;
				FatalError("Divergence should be close to 0, maybe there are some numerical round-off errors?");
			}

			/**
			 * TEST for correct vorticity
			 */
			/*
			 * Setup relative vorticity
			 */
			SphereData_Physical vortg(sphereDataConfig);
			vortg.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					// Eq. 94, Williamson TC paper

					// absolute vorticity, but we like the relative one
					//o_data = (2.0*u0/a + 2.0*omega)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));

					// relative vorticity
					o_data = (2.0*u0/a)*(-std::cos(lon)*std::cos(phi)*std::sin(alpha) + std::sin(phi)*std::cos(alpha));
				}
			);

			double vort_diff = (o_vort.getSphereDataPhysical() - vortg).physical_reduce_max_abs();
			if (vort_diff > 1e-12)
			{

				std::cout << "Vorticity difference: " << vort_diff << std::endl;
				FatalError("Vorticity fields differ (should be close to 0), maybe there are some numerical round-off errors?");
			}


			std::cout << "[MULE] geostrophic_balance_analytical_setup: 1" << std::endl;

			computeGeostrophicBalance_linear(
					o_vort,
					o_div,
					o_phi_pert
			);
		}
		else
		{
			FatalError("Benchmark not implemented");
		}
	}
};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SWESPHEREBENCHMARKSCOMBINED_HPP_ */
