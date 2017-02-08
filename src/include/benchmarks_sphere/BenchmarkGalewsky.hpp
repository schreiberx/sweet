/*
 * BenchmarkGalewsky.hpp
 *
 *  Created on: 16 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_BENCHMARKGALEWSKY_HPP_
#define SRC_BENCHMARKGALEWSKY_HPP_

#include <libmath/GaussQuadrature.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>

class BenchmarkGalewsky
{
	/**
	 * Stream jet data
	 */
	const double phi0 = M_PI/7.0;
	//const double phi0 = M_PI/30.0;
	const double phi1 = M_PI/2.0 - phi0;


	const double u_max = 80.0;
	const double h_avg = 10000.0;
	const double en = exp(-4.0/std::pow(phi1-phi0, 2.0));

	/**
	 * Bump data
	 */
	const double phi2 = M_PI/4.0;
	const double alpha = 1.0/3.0;
	const double beta = 1.0/15.0;
	const double h_hat = 120.0;


	SimulationVariables &simVars;

public:
	BenchmarkGalewsky(SimulationVariables &i_simVars)	:
		simVars(i_simVars)
	{

	}


private:
	static BenchmarkGalewsky** getPtrT()
	{
		// use this static pointer to allow using an existing quadrature code
		static BenchmarkGalewsky *t;
		return &t;
	}

	static BenchmarkGalewsky* T()
	{
		return *getPtrT();
	}

	double initial_condition_u(double lon, double phi)
	{
		if (phi <= phi0 || phi >= phi1)
			return 0;

		return u_max/en * exp(1.0/((phi-phi0)*(phi-phi1)));
	};

	void initial_condition_u(double lon, double lat, double &o_data)
	{
		o_data = initial_condition_u(lon, lat);
	};


	static double to_int_fun(double phi)
	{
		BenchmarkGalewsky* t = T();

		double u_phi = t->initial_condition_u(0, phi);

		return t->simVars.sim.earth_radius*u_phi*
				(2.0*t->simVars.sim.coriolis_omega*std::sin(phi)+(std::tan(phi)/t->simVars.sim.earth_radius)*u_phi);
	};

	double error_threshold = 1.e-13;

	double integrate_fun(
			double int_start,
			double int_end
	)
	{
#if 0
		static AdaptiveIntegrator<double(const double)> integration_helper;
		return integration_helper.integrate(to_int_fun, int_start, int_end, error_threshold);
#else
		return GaussQuadrature::integrate5_intervals_adaptive_linear<double>(int_start, int_end, to_int_fun, error_threshold);
#endif
	}


public:
	void setup_initial_h(SphereData &o_h)
	{
		assert(h_avg == 10000);

		// use this static pointer to allow using an existing quadrature code
		*getPtrT() = this;
		const SphereDataConfig *sphereDataConfig = o_h.sphereDataConfig;


		/*
		 * Initialization of U and V
		 */


		/*
		 * Initialization of H
		 *
		 * Metric correction terms based on John Thuburn's code
		 */
		const unsigned short nlat = sphereDataConfig->physical_num_lat;
		double *hg_cached =  new double[nlat];

		double h_area = 0;
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

			double hg = -integrate_fun(int_start, int_end);
			//hg = (int_end+int_start)*0.5;
			hg_cached[j] = hg;

			/*
			 * cos scaling is required for 2D sphere coverage at this latitude
			 *
			 * metric term which computes the area coverage of each point
			 */
			// use integrated average as below instead of the following formulation
			// double mterm = cos((int_start+int_end)*0.5);
			double mterm = (sin(int_end)-sin(int_start))*2.0*M_PI;
			assert(mterm > 0);

			hg_sum += hg*mterm;
			h_area += mterm;

			int_start = int_end;
		}
		j--;

		for (; j >= 0; j--)
		{
			double int_end = sphereDataConfig->lat[j];
			int_delta = int_end - int_start;
			assert(int_delta > 0);

			double hg = hg_cached[j+1] - integrate_fun(int_start, int_end);
			//hg = (int_end+int_start)*0.5;
			hg_cached[j] = hg;

			// metric term which computes the area coverage of each point
			double mterm = (sin(int_end)-sin(int_start))*2.0*M_PI;

			hg_sum += hg*mterm;
			h_area += mterm;

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
			double mterm = (sin(int_end)-sin(int_start))*2.0*M_PI;

			double hg = hg_cached[0] - integrate_fun(int_start, int_end);
			//hg = (int_end+int_start)*0.5;
			hg_sum += hg*mterm;
			h_area += mterm;
		}
		assert(h_area > 0);

		double h_sum = hg_sum / simVars.sim.gravitation;
		double h_comp_avg = h_sum / h_area;

		// shift to 10km
		for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
			hg_cached[j] = hg_cached[j]/simVars.sim.gravitation + (h_avg-h_comp_avg);

		// update data
		for (int i = 0; i < sphereDataConfig->physical_num_lon; i++)
			for (int j = 0; j < sphereDataConfig->physical_num_lat; j++)
				o_h.physical_space_data[i*sphereDataConfig->physical_num_lat + j] = hg_cached[j];

		o_h.physical_space_data_valid = true;
		o_h.spectral_space_data_valid = false;

		delete [] hg_cached;
	}


	void setup_initial_h_add_bump(SphereData &o_h)
	{
		o_h.physical_update_lambda(
				[&](double lambda, double phi, double &io_data)
				{
					io_data += h_hat*cos(phi)*exp(-pow((lambda-M_PI)/alpha, 2.0))*exp(-pow((phi2-phi)/beta, 2.0));
				}
		);
	}


	void setup_initial_u(SphereData &o_u)
	{
		o_u.physical_update_lambda(
				[&](double lon, double lat, double &o_data)
				{
					initial_condition_u(lon, lat, o_data);
				}
		);
	}


	void setup_initial_v(SphereData &o_v)
	{
		o_v.physical_set_zero();
	}

};


#endif /* SRC_BENCHMARKGALEWSKY_HPP_ */
