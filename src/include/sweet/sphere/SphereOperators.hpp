/*
 * SPHOperators.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHOPERATORS_HPP_
#define SPHOPERATORS_HPP_

#include <sweet/MemBlockAlloc.hpp>
#include "../sphere/SphereData.hpp"
#include "../sphere/SphereSPHIdentities.hpp"



class SphereOperators	:
		public SphereSPHIdentities
{
	friend SphereDataConfig;



public:
	/**
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
	SphereData diff_lon(
			const SphereData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SphereData out_sph_data(i_sph_data.sphereDataConfig);

		// compute d/dlambda in spectral space
#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int m = 0; m <= i_sph_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] = i_sph_data.spectral_space_data[idx]*std::complex<double>(0, m);
				idx++;
			}
		}
		out_sph_data.spectral_space_data_valid = true;
		out_sph_data.physical_space_data_valid = false;

		return out_sph_data;
	}



	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d mu f(lambda,mu)
	 *
	 * sqrt(1-mu*mu)*d/dmu P_n^m = ...
	 */
	SphereData diff_lat_mu(
			const SphereData &i_sph_data
	)	const
	{
		SphereData out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.physical_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					o_data /= (1.0-mu*mu);
				}
		);

		return out_sph_data;
	}



	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d phi f(lambda,mu)
	 */
	SphereData diff_lat_phi(
			const SphereData &i_sph_data
	)	const
	{
		return grad_lat(i_sph_data);
	}



	/**
	 * Compute gradient component along longitude (lambda)
	 */
	SphereData grad_lon(
			const SphereData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SphereData out_sph_data = diff_lon(i_sph_data);

		out_sph_data.request_data_physical();

		out_sph_data.physical_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					double cos_phi = sqrt(1.0-mu*mu);
					o_data /= cos_phi;
				}
		);

		return out_sph_data;
	}



	/**
	 * Divergence Operator along longitude for robert function formlation
	 *
	 * This computes
	 * 	1/cos^2(phi)  d/dlambda U
	 */
	SphereData robert_div_lon(
			const SphereData &i_sph_data
	)	const
	{
		// Entirely in spectral space
		SphereData out = diff_lon(i_sph_data);

		// Physical space
		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, double &o_data)
				{
					o_data /= cos_phi*cos_phi;
				}
			);

		return out;
	}



	/**
	 * Compute divergence along latitude for robert function formulation
	 *
	 * This computes
	 * 		d/dmu V
	 */
	SphereData robert_div_lat(
			const SphereData &i_sph_data
	)	const
	{
		/*
		 * Compute
		 *   cos^2(phi) * d/d mu  f(lambda,mu)
		 */
		// Entirely in spectral space
		SphereData out = spec_cosphi_squared_diff_lat_mu(i_sph_data);

		// Physical space
		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, double &o_data)
				{
					o_data /= cos_phi*cos_phi;
				}
			);

		return out;
	}



	/**
	 * Compute gradient component along longitude (lambda) for Robert function formulation
	 *
	 * This computes
	 * 		d/dlambda Phi
	 * with Phi the geopotential
	 */
	SphereData robert_grad_lon(
			const SphereData &i_sph_data
	)	const
	{
		// Entirely in spectral space
		return diff_lon(i_sph_data);
	}


	/**
	 * Compute gradient component along latitude for Robert function formulation
	 *
	 * This computes
	 * 		cos^2(phi) * d/dmu Phi
	 *
	 * with Phi the geopotential
	 */
	SphereData robert_grad_lat(
			const SphereData &i_sph_data
	)	const
	{
		// Entirely in spectral space
		//return spec_cosphi_squared_diff_lat_mu(i_sph_data);
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	inline
	SphereData spec_one_minus_sinphi_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)	const
	{
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	inline
	SphereData spec_cosphi_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)	const
	{
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	SphereData spec_one_minus_mu_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();
		SphereDataConfig *sphConfig = i_sph_data.sphereDataConfig;

		SphereData out_sph_data = SphereData(sphConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int m = 0; m <= i_sph_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] =	((-n+1.0)*R(n-1,m))*i_sph_data.spectral_get(n-1, m) +
												((n+2.0)*S(n+1,m))*i_sph_data.spectral_get(n+1, m);
				idx++;
			}
		}

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}



	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SphereData mu(
			const SphereData &i_sph_data
	)	const
	{
		SphereDataConfig *sphConfig = i_sph_data.sphereDataConfig;
		i_sph_data.request_data_spectral();

		SphereData out_sph_data = SphereData(sphConfig);


#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int m = 0; m <= i_sph_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] =
							R(n-1,m)*i_sph_data.spectral_get(n-1, m)
							+ S(n+1,m)*i_sph_data.spectral_get(n+1, m);

				idx++;
			}
		}

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}

	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SphereData mu2(
			const SphereData &i_sph_data
	)	const
	{
		SphereDataConfig *sphConfig = i_sph_data.sphereDataConfig;
		i_sph_data.request_data_spectral();

		SphereData out_sph_data = SphereData(sphConfig);


#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int m = 0; m <= i_sph_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] =
						+A(n-2,m)*i_sph_data.spectral_get(n-2, m)
						+B(n+0,m)*i_sph_data.spectral_get(n+0, m)
						+C(n+2,m)*i_sph_data.spectral_get(n+2, m)
						;
				idx++;
			}
		}

		out_sph_data.physical_space_data_valid = false;
		out_sph_data.spectral_space_data_valid = true;

		return out_sph_data;
	}


	/**
	 * Compute gradient component along latitude
	 */
	SphereData grad_lat(
			const SphereData &i_sph_data
	)	const
	{
		/*
		 * compute sin(theta)*d/d theta
		 * theta is the colatitude
		 *
		 * Hence, we have to
		 * 	first divide by sin(M_PI*0.5-phi) and
		 * 	second multiply by sqrt(1-mu*mu)
		 */
		SphereData out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.request_data_physical();
		out_sph_data.physical_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					//double phi = asin(mu);

					//o_data /= sin(M_PI*0.5-phi);
					//o_data /= ::cos(phi);
					o_data /= sqrt(1.0-mu*mu);
				}
		);

#if 0
		/**
		 * WARNING: Leave this code here
		 * We can see that the following operations would cancel out.
		 * Therefore this was commented.
		 */
		// undo the sin(theta) and multiply with sqrt(1-mu*mu)
		out_sph_data.request_data_physical();
		out_sph_data.physical_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					double phi = asin(mu);

					//o_data /= sin(M_PI*0.5-phi);
					o_data /= ::cos(phi);

					double cos_phi = std::sqrt((double)(1.0-mu*mu));
					o_data *= cos_phi;
				}
			);
#endif

		return out_sph_data;
	}



	/**
	 * Divergence Operator along longitude
	 *
	 * Identical to gradient operator along longitude
	 */
	SphereData div_lon(
			const SphereData &i_sph_data
	)	const
	{
		return grad_lon(i_sph_data);
	}



	/**
	 * Divergence Operator along latitude
	 *
	 * d(sqrt(1-mu*mu)*F)/dmu
	 */
	SphereData div_lat(
			const SphereData &i_sph_data
	)	const
	{
		SphereData out_sph_data(i_sph_data);

		// TODO: replace this with a recurrence identity
		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					//o_data *= cos(phi);
					o_data *= mu;
				}
			);

		// grad_lat = diff_lat_phi
		out_sph_data = grad_lat(out_sph_data);

		// undo the sin(theta) which is cos(phi)
		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data /= mu;
					//o_data /= cos(phi);
				}
			);

		return out_sph_data;
	}

#if 0
	/**
	 * Divergence Operator along latitude
	 *
	 * d(sqrt(1-mu*mu)*F)/dmu
	 */
	SphereData div_lat_TEST(
			const SphereData &i_sph_data
	)	const
	{
		SphereData out_sph_data(i_sph_data);

		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					//o_data *= cos(phi);
					o_data *= mu;
				}
			);

		// grad_lat = diff_lat_phi
		out_sph_data = grad_lat(out_sph_data);
#if 0
		// undo the sin(theta) which is cos(phi)
		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data /= mu;
					//o_data /= cos(phi);
				}
			);
#endif
		return out_sph_data;
	}
#endif



	/**
	 * Laplace operator
	 */
	SphereData laplace(
			const SphereData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SphereData out_sph_data(i_sph_data);

		out_sph_data.spectral_update_lambda(
				[](int n, int m, std::complex<double> &o_data)
				{
					o_data *= -(double)n*((double)n+1.0);
				}
			);

		return out_sph_data;
	}


public:
	/**
	 * Compute vorticity
	 *
	 * \eta = grad_lat(V_lon) - grad_lon(V_lat)
	 */
	SphereData vort(
			const SphereData &i_lon,
			const SphereData &i_lat
	)	const
	{
		return grad_lon(i_lat) - grad_lat(i_lon);
	}



public:
	/**
	 * Compute divergence
	 *
	 * \delta = div_lon(i_lon) + div_lan(i_lan)
	 */
	SphereData div(
			const SphereData &i_lon,
			const SphereData &i_lat
	)	const
	{
		return div_lon(i_lon) + div_lat(i_lat);
	}
};






#endif /* SPHOPERATORS_HPP_ */
