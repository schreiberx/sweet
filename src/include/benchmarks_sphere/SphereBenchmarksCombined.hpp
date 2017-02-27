/*
 * BenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>
#include <benchmarks_sphere/BenchmarkGalewsky.hpp>
#include <benchmarks_sphere/BenchmarkGaussianDam.hpp>

class SphereBenchmarksCombined
{
public:
	static
	void setupInitialConditions(
			SphereData &o_h,
			SphereData &o_u,
			SphereData &o_v,
			SimulationVariables &io_simVars,
			SphereOperators &i_op
	)
	{
		BenchmarkGalewsky benchmarkGalewsky(io_simVars);

		if (io_simVars.setup.benchmark_scenario_id < 0)
		{
			std::cout << std::endl;
			std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
			std::cout << "Available benchmark scenarios:" << std::endl;
			std::cout << "	0: Dummy" << std::endl;
			std::cout << "	2: Use Gaussian bump initial conditions (0, pi/3)" << std::endl;
			std::cout << "	3: Use Gaussian bump initial conditions (pi/3, pi/3)" << std::endl;
			std::cout << "	4: Use Gaussian bump initial conditions (-pi/2, pi/4)" << std::endl;
			std::cout << "	10: Use geostrophic balance test case" << std::endl;
			std::cout << "	11: Williamson test benchmark 1 (DIV formulation)" << std::endl;
			std::cout << "	12: Williamson test benchmark 1 (GRAD formulation)" << std::endl;
			std::cout << "	100: Galweski" << std::endl;
			std::cout << "	101: Galweski - geostrophic case including non-linear parts" << std::endl;
			std::cout << "	200: h=h0, u=0, v=0" << std::endl;
			std::cout << std::endl;
			FatalError("Benchmark scenario not selected");
		}

		if (io_simVars.setup.benchmark_scenario_id == 0)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 2)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, 0, M_PI/3.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 3)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 4)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, -M_PI, M_PI/4.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 5)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, -M_PI, M_PI/4.0, 100.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 6)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, M_PI/3.0, M_PI/3.0, 20.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 7)
		{
			BenchmarkGaussianDam::setup_initial_conditions_gaussian(o_h, o_u, o_v, io_simVars, -M_PI, M_PI/4.0, 100.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 8)
		{
			/*
			 * PAUL N. SWARZTRAUBER
			 * Shallow Water Flow on the Sphere
			 * 2004
			 */
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;

			if (io_simVars.sim.coriolis_omega != 0)
				io_simVars.sim.coriolis_omega = 7.292e-5;

			io_simVars.sim.gravitation = 9.80616;
			io_simVars.sim.earth_radius = 6.37122e6;
			io_simVars.sim.h0 = 29400.0;

			o_u.spectral_set_zero();
			o_v.spectral_set_zero();

			double a = io_simVars.sim.earth_radius;
			double A = 6000.0;
			double alpha = 10;

			o_h.physical_update_lambda(
				[&](double i_lambda, double i_phi, double &io_data)
				{
					double x = a*std::cos(i_phi)*std::cos(i_lambda);
					double y = a*std::cos(i_phi)*std::sin(i_lambda);
					double z = a*std::sin(i_phi);

					double d = std::sqrt(x*x+y*y+(z-a)*(z-a));

					io_data = io_simVars.sim.h0 + A*std::exp(-alpha*(d/a)*(d/a));
				}
			);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 10)
		{
			/*
			 * Williamson test case
			 */
			double inv_r = 1.0/io_simVars.sim.earth_radius;

			o_v.spectral_set_zero();

			o_h.physical_update_lambda(
				[&](double i_lon, double i_lat, double &io_data)
				{
					io_data = io_simVars.sim.earth_radius*io_simVars.sim.coriolis_omega*std::cos(i_lat)*std::cos(i_lat)/io_simVars.sim.gravitation;
				}
			);

			if (io_simVars.misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						io_data = std::cos(i_lat)*std::cos(i_lat);
					}
				);
			}
			else
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						//io_data = simVars.sim.earth_radius*2.0*simVars.sim.coriolis_omega*std::cos(i_lat)/simVars.sim.gravitation;
						io_data = std::cos(i_lat);
					}
				);
			}

			if (io_simVars.timecontrol.current_simulation_time == 0)
			{
				double h_max_error, u_max_error, v_max_error;

				if (io_simVars.misc.sphere_use_robert_functions)
				{
					h_max_error = (-inv_r*i_op.robert_div_lon(o_u) -inv_r*i_op.robert_div_lat(o_v)).physical_reduce_max_abs();
					u_max_error = (-inv_r*i_op.robert_grad_lon(o_h*io_simVars.sim.gravitation) + 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_v)).physical_reduce_max_abs();
					v_max_error = (-inv_r*i_op.robert_grad_lat(o_h*io_simVars.sim.gravitation) - 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_u)).physical_reduce_max_abs();
				}
				else
				{
					h_max_error = (-inv_r*i_op.div_lon(o_u) -inv_r*i_op.div_lat(o_v)).physical_reduce_max_abs();
					u_max_error = (-inv_r*i_op.grad_lon(o_h*io_simVars.sim.gravitation) + 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_v)).physical_reduce_max_abs();
					v_max_error = (-inv_r*i_op.grad_lat(o_h*io_simVars.sim.gravitation) - 2.0*io_simVars.sim.coriolis_omega*i_op.mu(o_u)).physical_reduce_max_abs();
				}
				std::cout << "h_max_error for geostrophic balance case: " << h_max_error << std::endl;
				std::cout << "u_max_error for geostrophic balance case: " << u_max_error << std::endl;
				std::cout << "v_max_error for geostrophic balance case: " << v_max_error << std::endl;
			}
		}
		else if (io_simVars.setup.benchmark_scenario_id == 11)
		{
			/*
			 * See Williamson test case, eq. (75), (76)
			 *
			 * phi_t = DIV (u phi)
			 */

			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;

			io_simVars.sim.coriolis_omega = 7.292e-5;
			io_simVars.sim.gravitation = 9.80616;
			io_simVars.sim.earth_radius = 6.37122e6;
			io_simVars.sim.h0 = 1000.0;

			double lambda_c = 3.0*M_PI/2.0;
			double theta_c = 0;
			double a = io_simVars.sim.earth_radius;

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
						io_data = io_simVars.sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
					else
						io_data = 0;
				}
			);


			if (io_simVars.misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
							u0*(
									std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
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
									std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
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
								std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
								std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
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
								std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
							);
					}
				);
			}

			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Storing advection in output velocities !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;


			io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);

			std::cout << "advection_rotation_angle: " << io_simVars.setup.advection_rotation_angle << std::endl;
		}
		else if (io_simVars.setup.benchmark_scenario_id == 12)
		{
			/**
			 * See Williamson test case, eq. (75), (76)
			 *
			 * phi_t = phi GRAD u
			 */

			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;

			io_simVars.sim.coriolis_omega = 7.292e-5;
			io_simVars.sim.gravitation = 9.80616;
			io_simVars.sim.earth_radius = 6.37122e6;
			io_simVars.sim.h0 = 1000.0;

			double lambda_c = 3.0*M_PI/2.0;
			double theta_c = 0;
			double a = io_simVars.sim.earth_radius;

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
						io_data = io_simVars.sim.h0/2.0*(1.0+std::cos(M_PI*r/R));
					else
						io_data = 0;
				}
			);


			if (io_simVars.misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda(
					[&](double i_lon, double i_lat, double &io_data)
					{
						double i_theta = i_lat;
						double i_lambda = i_lon;
						io_data =
								u0*(
									std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
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
									std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
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
									std::cos(i_theta)*std::cos(io_simVars.setup.advection_rotation_angle) +
									std::sin(i_theta)*std::cos(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
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
									std::sin(i_lambda)*std::sin(io_simVars.setup.advection_rotation_angle)
							);
					}
				);
			}

			o_h.physical_truncate();
			o_u.physical_truncate();
			o_v.physical_truncate();

			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Storing advection in output velocities !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;


			io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);

			std::cout << "advection_rotation_angle: " << io_simVars.setup.advection_rotation_angle << std::endl;
		}
		else if (io_simVars.setup.benchmark_scenario_id == 20)
		{
			o_h.physical_set_all_value(io_simVars.sim.h0);
			o_u.physical_set_all_value(0);
			o_v.physical_set_all_value(0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 1 || io_simVars.setup.benchmark_scenario_id == 100 || io_simVars.setup.benchmark_scenario_id == 101)
		{

#if 0

			benchmarkGalewsky.setup_initial_h(o_h);
#if 0
			std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
			std::cout << "TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!" << std::endl;
			/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
			/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
			/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
			/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
			/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
			/// TODOOOOOOOOOOOOOOOOOOOOOOO: REMOVE ME!!!!!!!!!!!!!!!!!!
			double foo[100] = {
					9055.493745200144986,  9055.465517385822750,  9055.482442174412427,  9055.485353917110842,  9055.470071095474850,  9055.482225177744112,  9055.483738383378295,  9055.470977648797998,  9055.482690202114100,  9055.483017123375248,  9055.470682326214956,  9055.483931642287644,  9055.482260546816178,  9055.469577987099910,  9055.486730172287935,  9055.479593312960787,  9055.470174185988071,  9055.506848781622466,  9056.095270992454971,  9061.214200969297963,  9081.885216784596196,  9134.617053737521928,  9233.437475196804371,  9380.388945759646958,  9560.630706589838155,  9746.223428225062889,  9907.553693800491601,  10025.289612307593416,  10095.839141933991414,  10128.926698702083740,  10140.037739517365480,  10142.276285688618373,  10142.446739097140380,  10142.457887177219163,  10142.460637987251175,  10142.454087103898928,  10142.457439067522500,  10142.458913848522570,  10142.455881685920758,  10142.456781481529106,  10142.458514041243689,  10142.456551124038015,  10142.456497322984433,  10142.458327286207350,  10142.456871132475499,  10142.456378846281950,  10142.458196108580523,  10142.457058206358852,  10142.456330244574929,  10142.458089350835508,  10142.457184941487867,  10142.456312120679286,  10142.457997722302025,  10142.457280603905019,  10142.456307374064636,  10142.457916977255081,  10142.457358891591866,  10142.456308088640071,  10142.457844461576315,  10142.457427072755308,  10142.456310368295817,  10142.457778217069063,  10142.457489475389593,  10142.456312155902197,  10142.457716672483002,  10142.457549000710060,  10142.456312239946783,  10142.457658493196504,  10142.457607867550905,  10142.456309759470969,  10142.457602479804336,  10142.457668042088699,  10142.456303918996127,  10142.457547479152709,  10142.457731553769918,  10142.456293774963342,  10142.457492286770503,  10142.457800817175666,  10142.456278003344778,  10142.457435515911129,  10142.457879081055580,  10142.456254545904812,  10142.457375387479260,  10142.457971214398640,  10142.456219940253504,  10142.457309336608887,  10142.458085327591107,  10142.456167834954613,  10142.457233155662834,  10142.458236684131407,  10142.456085108484331,  10142.457138761870738,  10142.458459248211511,  10142.455939082734403,  10142.457006667080350,  10142.458852741790906,  10142.455614339081876,  10142.456765204698968,  10142.459972075210317,  10142.453909706555351
			};

			o_h.physical_update_lambda_array(
					[&](int i_lon, int i_lat, double &io_data)
					{
						io_data = foo[i_lat];
					}
			);
#endif

			if (io_simVars.setup.benchmark_scenario_id != 101)
				benchmarkGalewsky.setup_initial_h_add_bump(o_h);

			benchmarkGalewsky.setup_initial_u(o_u);
			benchmarkGalewsky.setup_initial_v(o_v);
#if 0
			std::cout << o_u.sphereDataConfig->shtns->nlm << " " << o_u.sphereDataConfig->spectral_array_data_number_of_elements << std::endl;
			for (int i = 0; i < o_u.sphereDataConfig->physical_num_lat; i++)
				std::cout << i << ": " << o_u.physical_get(0, i) << std::endl;
			o_u.physical_file_write("o_u_0_physical.csv");
			o_u.request_data_spectral();
			o_u.physical_file_write("o_u_1_spectral_physical.csv");

			std::cout << "----" << std::endl;
			for (int i = 0; i < o_u.sphereDataConfig->physical_num_lat; i++)
				std::cout << i << ": " << o_u.physical_get(0, i) << std::endl;
#endif

			if (io_simVars.misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda_cosphi_grid(
					[&](double i_lon, double i_cosphi, double &io_data)
					{
						io_data *= i_cosphi;
					}
				);

				o_v.physical_update_lambda_cosphi_grid(
					[&](double i_lon, double i_cosphi, double &io_data)
					{
						io_data *= i_cosphi;
					}
				);
			}
/*
			o_h.spectral_truncate();
			o_u.spectral_truncate();
			o_v.spectral_truncate();
*/

#else

			/*
			 * Initialization code from Jeffrey Whitaker
			 */

			double umax = 80.;
			double phi0 = M_PI/7.;
			double phi1 = 0.5*M_PI - phi0;
			double phi2 = 0.25*M_PI;
			double en = std::exp(-4.0/std::pow((phi1-phi0), 2.0));
			double alpha = 1./3.;
			double beta = 1./15.;
			double hamp = 120.;

			SphereDataPhysical f(o_h.sphereDataConfig);
			f.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = 2.0*io_simVars.sim.coriolis_omega*mu;
				}
			);

			SphereDataPhysical vg(o_h.sphereDataConfig);
			vg.physical_set_zero();

		    //u1 = (umax/en)*np.exp(1./((x.lats-phi0)*(x.lats-phi1)))
			SphereDataPhysical u1(o_h.sphereDataConfig);
			u1.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					o_data = umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
				}
			);

			SphereDataPhysical ug = u1;
			ug.physical_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					if (phi >= phi1 || phi <= phi0)
						o_data = 0;
				}
			);

			SphereDataPhysical hbump(o_h.sphereDataConfig);

			if (io_simVars.setup.benchmark_scenario_id != 101)
			{
				hbump.physical_update_lambda(
					[&](double lon, double phi, double &o_data)
					{
						o_data = hamp*std::cos(phi)*std::exp(-(lon/alpha)*(lon/alpha))*std::exp(-(phi2-phi)*(phi2-phi)/beta);
					}
				);
			}
			else
			{
				hbump.physical_set_zero();
			}



			SphereData vrtspec(o_h.sphereDataConfig);
			SphereData divspec(o_h.sphereDataConfig);

			i_op.uv_to_vortdiv(ug, vg, vrtspec, divspec);


			SphereDataPhysical vrtg = vrtspec.getSphereDataPhysical();

			SphereDataPhysical tmpg1 = ug*(vrtg+f);
			SphereDataPhysical tmpg2 = vg*(vrtg+f);

			SphereData tmpspec1(o_h.sphereDataConfig);
			SphereData tmpspec2(o_h.sphereDataConfig);
			i_op.uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

			tmpspec2 = (0.5*(ug*ug+vg*vg));

			SphereData phispec = i_op.inv_laplace(tmpspec1) - tmpspec2;

			io_simVars.sim.h0 = 10000.0;
			SphereDataPhysical phig = io_simVars.sim.gravitation*(hbump+io_simVars.sim.h0) + phispec.getSphereDataPhysical();

			phispec = phig;
			phispec.spectral_truncate();

			////////////////////////////////

			i_op.vortdiv_to_uv(vrtspec, divspec, ug, vg);

			o_h = phispec/io_simVars.sim.gravitation;
			o_u = ug;
			o_v = vg;

			if (io_simVars.misc.sphere_use_robert_functions)
			{
				o_u.physical_update_lambda_cosphi_grid(
					[&](double lon, double phi, double &o_data)
					{
						o_data *= phi;
					}
				);

				o_v.physical_update_lambda_cosphi_grid(
					[&](double lon, double phi, double &o_data)
					{
						o_data *= phi;
					}
				);
			}
#endif


			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;

			/// Setup Galewski parameters
			io_simVars.sim.coriolis_omega = 7.292e-5;
			io_simVars.sim.gravitation = 9.80616;
			io_simVars.sim.earth_radius = 6.37122e6;
			io_simVars.sim.h0 = benchmarkGalewsky.h_avg;

			io_simVars.misc.output_time_scale = 1.0/(60.0*60.0);
		}
		else if (io_simVars.setup.benchmark_scenario_id == 200)
		{
			o_h.physical_set_all_value(io_simVars.sim.h0);
			o_u.spectral_set_zero();
			o_v.spectral_set_zero();
		}
		else
		{
			FatalError("Unknown scenario id");
		}
	}

};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
