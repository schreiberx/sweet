/*
 * Author: Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_NAIR_LAURITZEN_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_NAIR_LAURITZEN_HPP_

#include "BenchmarksSphereAdvection_interface.hpp"
#include <ostream>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereData_Config.hpp>

#include <sweet/SWEETVectorMath.hpp>


class BenchmarksSphereAdvection_nair_lauritzen_sl	: public BenchmarksSphereAdvection_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


public:
	BenchmarksSphereAdvection_nair_lauritzen_sl()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "nair_lauritzen_case_test"	||
				benchmark_name == "nair_lauritzen_case_1"	||
				benchmark_name == "nair_lauritzen_case_2"	||
				benchmark_name == "nair_lauritzen_case_2_k0.5"	||
				benchmark_name == "nair_lauritzen_case_3"	||
				benchmark_name == "nair_lauritzen_case_4"	||
				benchmark_name == "nair_lauritzen_case_4_ext_2"	||
				benchmark_name == "nair_lauritzen_case_4_ext_3" ||
				false
		;
	}


	void setup(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
	)
	{
		simVars = i_simVars;
		ops = i_ops;
	}


	bool has_time_varying_state()
	{
		return true;
	}


	std::string get_help()
	{
		std::ostringstream stream;
		stream << " * NAIR LAURITZEN SL TEST CASES:" << std::endl;

		stream << "    + 'nair_lauritzen_case_test'" << std::endl;
		stream << "    + 'nair_lauritzen_case_1'" << std::endl;
		stream << "    + 'nair_lauritzen_case_2'" << std::endl;
		stream << "    + 'nair_lauritzen_case_2_k0.5'" << std::endl;
		stream << "    + 'nair_lauritzen_case_3'" << std::endl;
		stream << "    + 'nair_lauritzen_case_4'" << std::endl;
		stream << "    + 'nair_lauritzen_case_4_ext_2'" << std::endl;
		stream << "    + 'nair_lauritzen_case_4_ext_3'" << std::endl;

		return stream.str();
	}


	void get_initial_state(
		std::vector<SphereData_Spectral*> &o_prognostic_fields,
		SphereData_Physical &o_u,
		SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 1, "Only scalar field supported for this benchmark!");

		get_initial_state(*o_prognostic_fields[0], o_u, o_v);
	}

	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Physical &o_u,
		SphereData_Physical &o_v
	)
	{
		/*
		 * Time-varying benchmark case 4 from
		 *
		 * R. Nair, P. Lauritzen "A class of deformational flow
		 * test cases for linear transport problems on the sphere"
		 */

		/*
		 * Setup parameters
		 */

		// use the radius from the command line parameter since this is should work flawless
		//simVars->sim.sphere_radius = 6.37122e6;
		simVars->sim.h0 = 1.0;				// h_max

		if (std::isinf(simVars->timecontrol.max_simulation_time))
			simVars->timecontrol.max_simulation_time = 12*24*60*60;		// default: 12 days

		// update operators
		ops->setup(ops->sphereDataConfig, &(simVars->sim));

		double i_lambda0 = M_PI/3;
		double i_theta0 = M_PI;
		double i_lambda1 = -M_PI/3;
		double i_theta1 = M_PI;

		if (benchmark_name == "nair_lauritzen_case_test")
		{
			i_lambda0 = M_PI;
			i_theta0 = M_PI/3;
			i_lambda1 = M_PI;
			i_theta1 = -M_PI/3;
		}
		else if (benchmark_name == "nair_lauritzen_case_1")
		{
			i_lambda0 = M_PI;
			i_theta0 = M_PI/3;
			i_lambda1 = M_PI;
			i_theta1 = -M_PI/3;
		}
		else if (
				benchmark_name == "nair_lauritzen_case_2" ||
				benchmark_name == "nair_lauritzen_case_4" ||
				benchmark_name == "nair_lauritzen_case_4_ext_2"
		)
		{
			i_lambda0 = 5*M_PI/6.0;
			i_theta0 = 0;
			i_lambda1 = 7*M_PI/6.0;
			i_theta1 = 0;
		}
		else if (
				benchmark_name == "nair_lauritzen_case_3" ||
				benchmark_name == "nair_lauritzen_case_4_ext_3"
		)
		{
			i_lambda0 = 3*M_PI/4.0;
			i_theta0 = 0;
			i_lambda1 = 5*M_PI/4.0;
			i_theta1 = 0;
		}


#if 0

		/**
		 * Initial condition for Cosine bell
		 *
		 * DO NOT USE THIS, SINCE IT'S NOT SUFFICIENTLY SMOOTH FOR THE CONVERGENCE BENCHMARKS!!!
		 *
		 * (Section 3.1.1)
		 */
		SphereData_Physical phi_pert_phys_1(sphereDataConfig);
		{
			// Cosine bells

			SphereData_Physical phi_pert_phys(sphereDataConfig);
			phi_pert_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// Constants
					double r = 0.5;
					double b = 0.1;
					double c = 0.9;
					double pi = M_PI;

					// eq between 12 and 13
					double r0 = std::acos(std::sin(i_theta0)*std::sin(i_theta) + std::cos(i_theta0)*(std::cos(i_theta)*std::cos(i_lambda-i_lambda0)));
					double r1 = std::acos(std::sin(i_theta1)*std::sin(i_theta) + std::cos(i_theta1)*(std::cos(i_theta)*std::cos(i_lambda-i_lambda1)));

					io_data = 0;
					// eq. 13
					if (r0 < r)
						io_data = b + c*simVars->sim.h0/2.0*(1.0 + std::cos(pi*r0/r));
					else if (r1 < r)
						io_data = b + c*simVars->sim.h0/2.0*(1.0 + std::cos(pi*r1/r));
					else
						io_data = b;
				}
			);

			o_phi_pert = phi_pert_phys;
		}

#else

		double b0 = 5;
		//double b0 = 20;

		/**
		 * Initial condition for smooth scalar field (Gaussian bump)
		 *
		 * (Section 3.1.2)
		 */
		SphereData_Physical phi_pert_phys_1(ops->sphereDataConfig);
		{
			// Bump 1

			// Caption Figure 1
			double x0[3];
			SWEETVectorMath::point_latlon_to_cartesian__scalar(i_lambda0, i_theta0, x0[0], x0[1], x0[2]);

			phi_pert_phys_1.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					double x[3];
					SWEETVectorMath::point_latlon_to_cartesian__scalar(i_lambda, i_theta, x[0], x[1], x[2]);

					double d =	(x[0] - x0[0])*(x[0] - x0[0]) +
								(x[1] - x0[1])*(x[1] - x0[1]) +
								(x[2] - x0[2])*(x[2] - x0[2]);

					io_data = std::exp(-b0*d);

					io_data *= simVars->sim.h0;
				}
			);
		}

		SphereData_Physical phi_pert_phys_2(ops->sphereDataConfig);
		{
			// Bump 2

			// Caption Figure 1
			double x0[3];
			SWEETVectorMath::point_latlon_to_cartesian__scalar(i_lambda1, i_theta1, x0[0], x0[1], x0[2]);

			phi_pert_phys_2.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
				double x[3];
				SWEETVectorMath::point_latlon_to_cartesian__scalar(i_lambda, i_theta, x[0], x[1], x[2]);

				double d =	(x[0] - x0[0])*(x[0] - x0[0]) +
							(x[1] - x0[1])*(x[1] - x0[1]) +
							(x[2] - x0[2])*(x[2] - x0[2]);

				io_data = std::exp(-b0*d);

				io_data *= simVars->sim.h0;
				}
			);
		}
/*

		if (o_phi_pert_phys != nullptr)
			*o_phi_pert_phys = phi_pert_phys_1 + phi_pert_phys_2;
*/

		o_phi_pert.loadSphereDataPhysical(phi_pert_phys_1 + phi_pert_phys_2);
#endif

		get_varying_velocities(o_u, o_v, 0);

		return;
	}


	/*
	 * Update fields for time-varying benchmarks
	 */
	void get_varying_velocities(
			SphereData_Physical &o_u_phys,
			SphereData_Physical &o_v_phys,
			double i_timestamp = 0
	)
	{
		/*********************************************************************
		 * Time-varying benchmark cases
		 *
		 * R. Nair, P. Lauritzen "A class of deformational flow
		 * test cases for linear transport problems on the sphere"
		 */
		if (benchmark_name == "nair_lauritzen_case_test")
		{
			// time for total deformation
			double T = simVars->timecontrol.max_simulation_time;

			// velocity
			double u0 = 2.0*M_PI*simVars->sim.sphere_radius/T;

			// we set k to 2.4 (p. 5)
			double k = 0.5;

			o_u_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = k * std::pow(std::sin(i_lambda/2), 2.0) * std::sin(2.0*i_theta) * std::cos(M_PI*i_timestamp/T);
				}
			);

			o_v_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = k/2.0 * std::sin(i_lambda) * std::cos(i_theta) * std::cos(M_PI*i_timestamp/T);
				}
			);

			o_u_phys *= u0;
			o_v_phys *= u0;
			return;
		}

		if (benchmark_name == "nair_lauritzen_case_1")
		{
			// time for total deformation
			double T = 1;
			double t = i_timestamp;
			double pi = M_PI;

			// we set k to 2.4 (p. 5)
			double k = 2.4;

			/*
			 * Non-dimensionalize
			 */
			// t \in [0;1]
			t /= simVars->timecontrol.max_simulation_time;

			// velocity scalar for effects across entire simulation time (e.g. full revelation around sphere)
			// reference solution is computed with T = 5.0 and we resemble it here
			// using radius 1, T=5 would result in u0=1
			double u0 = simVars->sim.sphere_radius*5.0/simVars->timecontrol.max_simulation_time;

			using namespace ScalarDataArray_ops;

			o_u_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = k * pow2(sin(i_lambda/2.0)) * sin(2.0*i_theta) * cos(pi*t/T);
					io_data *= u0;
				}
			);

			o_v_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = k/2.0 * sin(i_lambda) * cos(i_theta) * cos(pi*t/T);
					io_data *= u0;
				}
			);

			return;
		}

		if (benchmark_name == "nair_lauritzen_case_2" || benchmark_name == "nair_lauritzen_case_2_k0.5")
		{
			// time for total deformation
			double T = 1;
			double t = i_timestamp;
			double pi = M_PI;

			// we set k to 2 (p. 5)
			double k = 2;

			/*
			 * Non-dimensionalize
			 */
			// t \in [0;1]
			t /= simVars->timecontrol.max_simulation_time;

			// velocity scalar for effects across entire simulation time (e.g. full revelation around sphere)
			// reference solution is computed with T = 5.0 and we resemble it here
			double u0 = simVars->sim.sphere_radius*5.0/simVars->timecontrol.max_simulation_time;

			using namespace ScalarDataArray_ops;

			if (benchmark_name == "nair_lauritzen_case_2_k0.5")
				k = 0.5;

			o_u_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = k * pow2(sin(i_lambda)) * sin(2.0*i_theta) * cos(pi*t/T);
					io_data *= u0;
				}
			);

			o_v_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = k * sin(2.0*i_lambda) * cos(i_theta) * cos(pi*t/T);
					io_data *= u0;
				}
			);

			return;
		}

		if (benchmark_name == "nair_lauritzen_case_3")
		{
			/*
			 * This is a divergent flow!!!
			 */

			// time for total deformation
			double T = 1;
			double t = i_timestamp;
			double pi = M_PI;

			// we set k to 2 (p. 7, "other parameters exactly as given in Case-2")
			double k = 2;

			/*
			 * Non-dimensionalize
			 */
			// t \in [0;1]
			t /= simVars->timecontrol.max_simulation_time;

			// velocity scalar for effects across entire simulation time (e.g. full revelation around sphere)
			// reference solution is computed with T = 5.0 and we resemble it here
			double u0 = simVars->sim.sphere_radius*5.0/simVars->timecontrol.max_simulation_time;

			using namespace ScalarDataArray_ops;

			o_u_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = -k * pow2(sin(i_lambda/2.0)) * sin(2.0*i_theta) * pow2(cos(i_theta)) * std::cos(pi*t/T);
					io_data *= u0;
				}
			);

			o_v_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					// time varying flow
					io_data = k/2.0 * sin(i_lambda) * pow3(cos(i_theta)) * cos(pi*t/T);
					io_data *= u0;
				}
			);
			return;
		}

		if (
			benchmark_name == "nair_lauritzen_case_4"	||
			benchmark_name == "nair_lauritzen_case_4_ext_2"
		)
		{

			// time for total deformation
			double T = 1;
			double t = i_timestamp;
			double pi = M_PI;

			// we set k to 2 (p. 5)
			double k = 2;

			/*
			 * Non-dimensionalize
			 */
			// t \in [0;1]
			t /= simVars->timecontrol.max_simulation_time;

			// velocity scalar for effects across entire simulation time (e.g. full revelation around sphere)
			// reference solution is computed with T = 5.0 and we resemble it here
			double u0 = simVars->sim.sphere_radius*5.0/simVars->timecontrol.max_simulation_time;

			double u0_rot = simVars->sim.sphere_radius/simVars->timecontrol.max_simulation_time;

			using namespace ScalarDataArray_ops;

			if (benchmark_name == "nair_lauritzen_case_2_k0.5")
				k = 0.5;

			o_u_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					io_data = 0;

					// washing machine
					double lambda_prime = i_lambda - 2.0*pi*t / T;
					io_data += u0 * k * pow2(sin(lambda_prime)) * sin(2.0*i_theta) * cos(pi*t/T);

					// add a constant zonal flow
					// non-dimensional version
					io_data += u0_rot*2.0*pi*cos(i_theta)/T;

				}
			);

			o_v_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					io_data = 0;

					// washing machine
					double lambda_prime = i_lambda - 2.0*pi*t / T;
					io_data += k * sin(2.0*lambda_prime) * cos(i_theta) * cos(pi*t/T);
					io_data *= u0;
				}
			);

			return;
		}

		if (benchmark_name == "nair_lauritzen_case_4_ext_3")
		{
			/*
			 * This is a modified case4 version to include the nonlinear divergence
			 */

			/*
			 * This is a divergent flow!!!
			 */

			// time for total deformation
			double T = 1;
			double t = i_timestamp;
			double pi = M_PI;

			// we set k to 2 (p. 7, "other parameters exactly as given in Case-2")
			double k = 2;

			/*
			 * Scale things up for the sphere
			 */
			// total simulation time
			T *= simVars->timecontrol.max_simulation_time;

			// velocity scalar for effects across entire simulation time (e.g. full revelation around sphere)
			// reference solution is computed with T = 5.0 and we resemble it here
			double u0 = simVars->sim.sphere_radius*5.0/T;

			using namespace ScalarDataArray_ops;

			o_u_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					io_data = 0;

					// washing machine
					double lambda_prime = i_lambda - 2.0*pi*t/T;
					io_data += -k * pow2(sin(lambda_prime/2.0)) * sin(2.0*i_theta) * pow2(cos(i_theta)) * std::cos(pi*t/T);
					io_data *= u0;

					// add a constant zonal flow
					// non-dimensional version
					//io_data += 2.0*pi*cos(i_theta)/T;
					io_data += 2.0*pi*cos(i_theta)*simVars->sim.sphere_radius/T;
				}
			);

			o_v_phys.physical_update_lambda(
				[&](double i_lambda, double i_theta, double &io_data)
				{
					io_data = 0;

					// washing machine
					double lambda_prime = i_lambda - 2.0*pi*t / T;
					io_data += k/2.0 * sin(lambda_prime) * pow3(cos(i_theta)) * cos(pi*t/T);
					io_data *= u0;
				}
			);

			return;
		}
	}



	/*
	 * Compute the departure points
	 */
	void sl_compute_departure_3rd_order(
			const ScalarDataArray &i_pos_lon_A,	///< longitude coordinate to compute the velocity for
			const ScalarDataArray &i_pos_lat_A,	///< latitude coordinate to compute the velocity for
			ScalarDataArray &o_pos_lon_D,		///< velocity along longitude
			ScalarDataArray &o_pos_lat_D,		///< velocity along latitude
			double i_dt,
			double i_timestamp_arrival			///< timestamp at arrival point
	)
	{
		if (benchmark_name == "nair_lauritzen_case_1")
		{
			/*********************************************************************
			 * R. Nair, P. Lauritzen "A class of deformational flow
			 * test cases for linear transport problems on the sphere"
			 *
			 * Page 8
			 */

			// time for total deformation
			double pi = M_PI;

			// we set k to 2.4 (p. 5)
			double k = 2.4;


			/*
			 * Get non-dimensionalize variables to stick to the benchmark
			 * whatever resolution we're using
			 *
			 * We rescale everything as if we would execute it with SWEET runtime parameters
			 *   -a 1
			 *   -t 5
			 */

			// Set full time interval to 5 (used in paper)
			double T = 5.0;

			// rescale time interval to [0;5] range
			double t = i_timestamp_arrival*5.0/simVars->timecontrol.max_simulation_time;

			// rescale time step size
			double dt = i_dt*5.0/simVars->timecontrol.max_simulation_time;


			double omega = pi/T;


			// To directly use notation from paper
			const ScalarDataArray &lambda = i_pos_lon_A;
			const ScalarDataArray &theta = i_pos_lat_A;


			// use convenient cos/sin/etc functions
			using namespace ScalarDataArray_ops;

			// eq. (37)
			ScalarDataArray u_tilde = 2.0*k*pow2(sin(lambda/2))*sin(theta)*cos(pi*t/T);

			// between eq. (34) and (35)
			ScalarDataArray v = k/2.0*sin(lambda)*cos(theta)*cos(omega*t);

			// eq. (35)
			o_pos_lon_D =
					lambda
					- dt*u_tilde
					- dt*dt*k*sin(lambda/2.0)*(
							sin(lambda/2)*sin(theta)*sin(omega*t)*omega
							- u_tilde*sin(theta)*cos(omega*t)*cos(lambda/2)
							- v*sin(lambda/2)*cos(theta)*cos(omega*t)
							);

			// eq. (36)
			o_pos_lat_D =
					theta
					- dt*v
					- dt*dt/4.0*k*(
							sin(lambda)*cos(theta)*sin(omega*t)*omega
							- u_tilde*cos(lambda)*cos(theta)*cos(omega*t)
							+ v*sin(lambda)*sin(theta)*cos(omega*t)
						);

			return;
		}

		SWEETError("TODO: Implement it for this test case");
	}


};

#endif
