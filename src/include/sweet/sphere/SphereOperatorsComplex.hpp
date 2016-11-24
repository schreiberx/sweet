/*
 * SPHOperatorsComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHEREOPERATORS_COMPLEX_HPP_
#define SPHEREOPERATORS_COMPLEX_HPP_

#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/SphereSPHIdentities.hpp>

class SphereOperatorsComplex	:
		public SphereSPHIdentities
{
	friend SphereDataConfig;

	/**
	 * Dummy constructor
	 */
public:
	SphereOperatorsComplex(
			SphereDataConfig *i_sphereDataConfig
	)
	{
	}


	/**
	 * Dummy constructor
	 */
public:
	SphereOperatorsComplex()
	{
	}



	/**
	 * Dummy setup
	 */
public:
	void setup(
			SphereDataConfig *i_sphereDataConfig,
			const double i_domain_size[2]
	)
	{
	}


public:
	/**
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
	static
	SphereDataComplex diff_lon(
			const SphereDataComplex &i_sph_data
	)
	{
		i_sph_data.request_data_spectral();

		SphereDataComplex out_sph_data(i_sph_data.sphereDataConfig);

		// compute d/dlambda in spectral space
#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
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
	static
	SphereDataComplex diff_lat_mu(
			const SphereDataComplex &i_sph_data
	)
	{

		SphereDataComplex out_sph_data = spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
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
	static
	SphereDataComplex diff_lat_phi(
			const SphereDataComplex &i_sph_data
	)
	{
		return grad_lat(i_sph_data);
	}



	/**
	 * Compute gradient component along longitude (lambda)
	 */
	static
	SphereDataComplex grad_lon(
			const SphereDataComplex &i_sph_data
	)
	{
		SphereDataComplex out_sph_data = diff_lon(i_sph_data);

		out_sph_data.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					double cos_phi = std::sqrt(1.0-mu*mu);
					o_data /= cos_phi;
				}
		);

		return out_sph_data;
	}


	static
	SphereDataComplex spectral_one_minus_mu_squared_diff_lat_mu(
			const SphereDataComplex &i_sph_data
	)
	{
		i_sph_data.request_data_spectral();

		const SphereDataConfig *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereDataComplex out_sph_data(sphereDataConfig);


#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				/**
				 * TODO: Optimize me!
				 */
				out_sph_data.spectral_space_data[idx] =
						((-n+1.0)*R(n-1,m))*i_sph_data.spectral_get(n-1, m) +
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
	static
	SphereDataComplex mu(
			const SphereDataComplex &i_sph_data
	)
	{
		const SphereDataConfig *sphConfig = i_sph_data.sphereDataConfig;
		i_sph_data.request_data_spectral();

		SphereDataComplex out_sph_data = SphereDataComplex(sphConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
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
	static
	SphereDataComplex mu2(
			const SphereDataComplex &i_sph_data
	)
	{
		const SphereDataConfig *sphConfig = i_sph_data.sphereDataConfig;
		i_sph_data.request_data_spectral();

		SphereDataComplex out_sph_data = SphereDataComplex(sphConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
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
	static
	SphereDataComplex grad_lat(
			const SphereDataComplex &i_sph_data
	)
	{
		/*
		 * compute sin(theta)*d/d theta
		 * theta is the colatitude
		 *
		 * Hence, we have to
		 * 	first divide by sin(M_PI*0.5-phi) and
		 * 	second multiply by sqrt(1-mu*mu)
		 */
		SphereDataComplex out_sph_data = spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
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
	static
	SphereDataComplex div_lon(
			const SphereDataComplex &i_sph_data
	)
	{
		return grad_lon(i_sph_data);
	}




	/**
	 * Divergence Operator along latitude
	 *
	 * d(sqrt(1-mu*mu)*F)/dmu
	 */
	static
	SphereDataComplex div_lat(
			const SphereDataComplex &i_sph_data
	)
	{
		SphereDataComplex out_sph_data(i_sph_data);

#if 1
		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					//o_data *= cos(phi);
					o_data *= mu;
				}
			);
#else
		out_sph_data.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					//o_data *= cos(phi);
					o_data *= std::sqrt(1.0-mu*mu);
				}
			);
#endif
		// grad_lat = diff_lat_phi

#if 1

		out_sph_data = spectral_one_minus_mu_squared_diff_lat_mu(out_sph_data);

		out_sph_data.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					//o_data /= mu;
					o_data /= (1.0-mu*mu);
				}
			);

#else
		out_sph_data = grad_lat(out_sph_data);

		// undo the sin(theta) which is cos(phi)
#if 0
		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					o_data /= mu;
					//o_data /= cos(phi);
				}
			);
#else
		out_sph_data.physical_update_lambda(
				[](double lambda, double phi, std::complex<double> &o_data)
				{
					//o_data /= mu;
					o_data /= cos(phi);
				}
			);
#endif
#endif

		return out_sph_data;
	}



	/**
	 * Divergence Operator along longitude for robert function formlation
	 *
	 * This computes
	 * 	1/cos^2(phi)  d/dlambda U
	 */
	static
	SphereDataComplex robert_div_lon(
			const SphereDataComplex &i_sph_data
	)
	{
		SphereDataComplex out = diff_lon(i_sph_data);

		// TODO: replace by colatitude function
		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, std::complex<double> &o_data)
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
	 *
	 * There's no other metric term involved!
	 */
	static
	SphereDataComplex robert_div_lat(
			const SphereDataComplex &i_sph_data
	)
	{
		/*
		 * Compute
		 *   cos^2(phi) * d/d mu  f(lambda,mu)
		 */
		SphereDataComplex out = spectral_cosphi_squared_diff_lat_mu(i_sph_data);

		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, std::complex<double> &o_data)
				{
					o_data /= cos_phi*cos_phi;
				}
			);

		return out;
	}



	static
	SphereDataComplex robert_cos2phi_div_lat(
			const SphereDataComplex &i_sph_data
	)
	{
		/*
		 * Compute
		 *   d/d mu  f(lambda,mu)
		 */
		SphereDataComplex out = spectral_cosphi_squared_diff_lat_mu(i_sph_data);

		return out;
	}



	/**
	 * Compute gradient component along longitude (lambda) for Robert function formulation
	 *
	 * This computes
	 * 		d/dlambda Phi
	 */
	static
	SphereDataComplex robert_grad_lon(
			const SphereDataComplex &i_sph_data
	)
	{
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
	static
	SphereDataComplex robert_grad_lat(
			const SphereDataComplex &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}


	inline
	static
	SphereDataComplex spectral_one_minus_sinphi_squared_diff_lat_mu(
			const SphereDataComplex &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	inline
	static
	SphereDataComplex spectral_cosphi_squared_diff_lat_mu(
			const SphereDataComplex &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * Laplace operator
	 */
	static
	SphereDataComplex laplace(
			const SphereDataComplex &i_sph_data
	)
	{
		i_sph_data.request_data_spectral();

		SphereDataComplex out(i_sph_data);

		out.spectral_update_lambda(
				[](int n, int m, std::complex<double> &o_data)
				{
					o_data *= -(double)n*((double)n+1.0);
				}
			);

		return out;
	}


public:
	/**
	 * Compute vorticity
	 *
	 * \eta = div_lat(V_lon) - div_lon(V_lat)
	 */
	static
	SphereDataComplex vort(
			const SphereDataComplex &i_lon,
			const SphereDataComplex &i_lat
	)
	{
		return div_lon(i_lat) - div_lat(i_lon);
	}


public:
	/**
	 * Compute vorticity
	 *
	 * \eta = grad_lat(V_lon) - grad_lon(V_lat)
	 */
	static
	SphereDataComplex robert_vort(
			const SphereDataComplex &i_lon,
			const SphereDataComplex &i_lat
	)
	{
		return robert_div_lon(i_lat) - robert_div_lat(i_lon);
	}



public:
	/**
	 * Compute divergence
	 *
	 * \delta = div_lon(i_lon) + div_lan(i_lan)
	 */
	static
	SphereDataComplex div(
			const SphereDataComplex &i_lon,
			const SphereDataComplex &i_lat
	)
	{
		return div_lon(i_lon) + div_lat(i_lat);
	}



public:
	/**
	 * Compute divergence
	 *
	 * \delta = div_lon(i_lon) + div_lan(i_lan)
	 */
	static
	SphereDataComplex robert_div(
			const SphereDataComplex &i_lon,
			const SphereDataComplex &i_lat
	)
	{
#if 1
		return robert_div_lon(i_lon) + robert_div_lat(i_lat);
#else
		// Only compute division once
		SphereDataComplex out = diff_lon(i_lon) + spectral_cosphi_squared_diff_lat_mu(i_lat);

		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, std::complex<double> &o_data)
				{
					o_data /= cos_phi*cos_phi;
				}
			);

		return out;
#endif
	}
};






#endif /* SPHOPERATORS_HPP_ */
