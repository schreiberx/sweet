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


/*
 * 0: Compute on-the-fly and element wise (independent from previous iteration)
 * 1: Start by direct computation and then use deltas [TODO: Not yet implemented]
 * 2: Use memory cached version
 */
#define SWEET_SPH_ON_THE_FLY_MODE	0


class SphereOperators	:
	public SphereSPHIdentities
{
	friend SphereDataConfig;

	SphereDataConfig *sphereDataConfig;

#if SWEET_SPH_ON_THE_FLY_MODE == 2
	std::vector<double> spec_one_minus_mu_squared_diff_lat_mu__1;
	std::vector<double> spec_one_minus_mu_squared_diff_lat_mu__2;

	std::vector<double> mu__1;
	std::vector<double> mu__2;

	std::vector<double> mu2__1;
	std::vector<double> mu2__2;
	std::vector<double> mu2__3;
#endif

	SphereData sqrt_one_minus_mu2;
	SphereData inv_one_minus_mu2;
	SphereData inv_sqrt_one_minus_mu2;


public:
	SphereOperators(
			SphereDataConfig *i_sphereDataConfig
	)	:
		sphereDataConfig(i_sphereDataConfig),
		sqrt_one_minus_mu2(i_sphereDataConfig),
		inv_one_minus_mu2(i_sphereDataConfig),
		inv_sqrt_one_minus_mu2(i_sphereDataConfig)
	{
		sqrt_one_minus_mu2.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data = std::sqrt(1.0-mu*mu);
				}
			);
		inv_one_minus_mu2.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data = 1.0/(1.0-mu*mu);
				}
			);
		inv_sqrt_one_minus_mu2.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data = 1.0/std::sqrt(1.0-mu*mu);
				}
			);


#if SWEET_SPH_ON_THE_FLY_MODE == 2
		std::size_t storage_size = sphereDataConfig->spectral_complex_array_data_number_of_elements;
		spec_one_minus_mu_squared_diff_lat_mu__1.resize(storage_size);
		spec_one_minus_mu_squared_diff_lat_mu__2.resize(storage_size);

		mu__1.resize(storage_size);
		mu__2.resize(storage_size);

		mu2__1.resize(storage_size);
		mu2__2.resize(storage_size);
		mu2__3.resize(storage_size);

		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				spec_one_minus_mu_squared_diff_lat_mu__1[idx] = (-n+1.0)*R(n-1,m);
				spec_one_minus_mu_squared_diff_lat_mu__2[idx] = (n+2.0)*S(n+1,m);

				mu__1[idx] = R(n-1,m);
				mu__2[idx] = S(n+1,m);

				mu2__1[idx] = A(n-2,m);
				mu2__2[idx] = B(n+0,m);
				mu2__3[idx] = C(n+2,m);

				idx++;
			}
		}
#endif
	}



public:
	/**
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
	SphereData diff_lon(
			const SphereData &i_sph_data
	)
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
	)
	{
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data)*inv_one_minus_mu2;
	}



	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d phi f(lambda,mu)
	 */
	SphereData diff_lat_phi(
			const SphereData &i_sph_data
	)
	{
		return grad_lat(i_sph_data);
	}



	/**
	 * Compute gradient component along longitude (lambda)
	 */
	SphereData grad_lon(
			const SphereData &i_sph_data
	)
	{
		i_sph_data.request_data_spectral();

		SphereData out_sph_data = diff_lon(i_sph_data);

		// physical space already requested if spectral space data is valid
		out_sph_data.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data /= std::sqrt(1.0-mu*mu);
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
	)
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
	)
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
	)
	{
		// Entirely in spectral space
		return diff_lon(i_sph_data);
	}



	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 */
	SphereData robert_grad_lon_M(
			const SphereData &i_sph_data
	)
	{
		SphereData retval = robert_grad_lon(i_sph_data);

		retval.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					double cos = std::cos(i_lat);
					io_data /= cos*cos;
				}
		);

		return retval;
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
	)
	{
		// Entirely in spectral space
		//return spec_cosphi_squared_diff_lat_mu(i_sph_data);
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 */
	SphereData robert_grad_lat_M(
			const SphereData &i_sph_data
	)
	{
		SphereData retval = robert_grad_lat(i_sph_data);

		retval.physical_update_lambda(
				[](double i_lon, double i_lat, double &io_data)
				{
					double cos = std::cos(i_lat);
					io_data /= cos*cos;
				}
		);

		return retval;
	}



	SphereData spec_one_minus_sinphi_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)
	{
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	SphereData spec_cosphi_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)
	{
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	SphereData spec_one_minus_mu_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)
	{
		i_sph_data.request_data_spectral();
		const SphereDataConfig *sphConfig = i_sph_data.sphereDataConfig;

		SphereData out_sph_data = SphereData(sphConfig);

#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int m = 0; m <= i_sph_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
#if SWEET_SPH_ON_THE_FLY_MODE == 0
				out_sph_data.spectral_space_data[idx] =	((-n+1.0)*R(n-1,m))*i_sph_data.spectral_get(n-1, m) +
												((n+2.0)*S(n+1,m))*i_sph_data.spectral_get(n+1, m);
#elif SWEET_SPH_ON_THE_FLY_MODE == 2
				out_sph_data.spectral_space_data[idx] =
						spec_one_minus_mu_squared_diff_lat_mu__1[idx]*i_sph_data.spectral_get(n-1, m)
						+ spec_one_minus_mu_squared_diff_lat_mu__2[idx]*i_sph_data.spectral_get(n+1, m);
#else
#	error "unsupported"
#endif
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
			const SphereData &i_sphere_data
	)
	{
		const SphereDataConfig *sphereDataConfig = i_sphere_data.sphereDataConfig;
		i_sphere_data.request_data_spectral();

		SphereData out_sph_data = SphereData(sphereDataConfig);


#if SWEET_THREADING
#pragma omp parallel for
#endif
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
#if SWEET_SPH_ON_THE_FLY_MODE == 0
				out_sph_data.spectral_space_data[idx] =
							R(n-1,m)*i_sphere_data.spectral_get(n-1, m)
							+ S(n+1,m)*i_sphere_data.spectral_get(n+1, m);
#elif SWEET_SPH_ON_THE_FLY_MODE == 2
				out_sph_data.spectral_space_data[idx] =
							mu__1[idx]*i_sphere_data.spectral_get(n-1, m)
							+ mu__2[idx]*i_sphere_data.spectral_get(n+1, m);
#else
#	error "unsupported"
#endif
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
	)
	{
		const SphereDataConfig *sphConfig = i_sph_data.sphereDataConfig;
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
#if SWEET_SPH_ON_THE_FLY_MODE == 0
				out_sph_data.spectral_space_data[idx] =
						+A(n-2,m)*i_sph_data.spectral_get(n-2, m)
						+B(n+0,m)*i_sph_data.spectral_get(n+0, m)
						+C(n+2,m)*i_sph_data.spectral_get(n+2, m)
						;
#elif SWEET_SPH_ON_THE_FLY_MODE == 2
				out_sph_data.spectral_space_data[idx] =
						+mu2__1[idx]*i_sph_data.spectral_get(n-2, m)
						+mu2__2[idx]*i_sph_data.spectral_get(n+0, m)
						+mu2__3[idx]*i_sph_data.spectral_get(n+2, m)
						;
#else
#	error "unsupported"
#endif
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
//	static
	SphereData grad_lat(
			const SphereData &i_sph_data
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
#if 0
		return spec_one_minus_mu_squared_diff_lat_mu(i_sph_data)*inv_sqrt_one_minus_mu2;
#else
		SphereData out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.request_data_physical();
		out_sph_data.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					//double phi = asin(mu);

					//o_data /= sin(M_PI*0.5-phi);
					//o_data /= ::cos(phi);
					o_data /= std::sqrt(1.0-mu*mu);
				}
		);

		return out_sph_data;
#endif

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
		return out_sph_data;
#endif

	}



	/**
	 * Divergence Operator along longitude
	 *
	 * Identical to gradient operator along longitude
	 */
//	static
	SphereData div_lon(
			const SphereData &i_sph_data
	)
	{
		return grad_lon(i_sph_data);
	}



	/**
	 * Divergence Operator along latitude
	 *
	 * d(sqrt(1-mu*mu)*F)/dmu
	 */
//	static
	SphereData div_lat(
			const SphereData &i_sph_data
	)
	{
#if 0
		SphereDataConfig sphereDataConfigExt;

		sphereDataConfigExt.setupAdditionalModes(i_sph_data.sphereDataConfig, 32, 32);

		SphereData data = i_sph_data.spectral_returnWithDifferentModes(&sphereDataConfigExt);

		SphereData sqrt_one_minus_mu2(&sphereDataConfigExt);
		sqrt_one_minus_mu2.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data = std::sqrt(1.0-mu*mu);
				}
			);

		SphereData inv_one_minus_mu2(&sphereDataConfigExt);
		inv_one_minus_mu2.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data = 1.0/(1.0-mu*mu);
				}
			);


		data = data*sqrt_one_minus_mu2;

		data = spec_one_minus_mu_squared_diff_lat_mu(data);

		data = data*inv_one_minus_mu2;

		return data.spectral_returnWithDifferentModes(i_sph_data.sphereDataConfig);

#else
		SphereData out_sph_data(i_sph_data);


#if 1
		out_sph_data = out_sph_data*sqrt_one_minus_mu2;

		out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(out_sph_data);

		out_sph_data = out_sph_data*inv_one_minus_mu2;

#else

		// TODO: replace this with a recurrence identity if possible
		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double comu, double &o_data)
				{
					//o_data *= cos(phi);
					o_data *= comu;
				}
			);

		// grad_lat = diff_lat_phi
		out_sph_data = grad_lat(out_sph_data);

		// undo the sin(theta) which is cos(phi)
		out_sph_data.physical_update_lambda_cogaussian_grid(
				[](double lambda, double comu, double &o_data)
				{
					o_data /= comu;
					//o_data /= cos(phi);
				}
			);

#endif

		return out_sph_data;
#endif
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
//	static
	SphereData laplace(
			const SphereData &i_sph_data
	)
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
	)
	{
		return div_lon(i_lat) - div_lat(i_lon);
	}


public:
	/**
	 * Compute vorticity
	 *
	 * \eta = robert_grad_lat(V_lon) - robert_grad_lon(V_lat)
	 */
	SphereData robert_vort(
			const SphereData &i_lon,
			const SphereData &i_lat
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
	SphereData div(
			const SphereData &i_lon,
			const SphereData &i_lat
	)
	{
		return div_lon(i_lon) + div_lat(i_lat);
	}


public:
	/**
	 * Compute divergence
	 *
	 * \delta = robert_div_lon(i_lon) + robert_div_lan(i_lan)
	 */
	SphereData robert_div(
			const SphereData &i_lon,
			const SphereData &i_lat
	)
	{
		return robert_div_lon(i_lon) + robert_div_lat(i_lat);
	}
};



#endif /* SPHOPERATORS_HPP_ */
