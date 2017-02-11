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
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>

/*
 * 0: Compute on-the-fly and element wise (independent from previous iteration)
 * 1: Start by direct computation and then use deltas [TODO: Not yet implemented]
 * 2: Use memory cached version
 */
#define SWEET_SPH_ON_THE_FLY_MODE	2


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


	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolver_inv_one_minus_mu2;


public:
	SphereOperators(
		SphereDataConfig *i_sphereDataConfig
	)
	{
		setup(i_sphereDataConfig);
	}


public:
	SphereOperators()
	{
	}



public:
	void setup(
		SphereDataConfig *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		sphSolver_inv_one_minus_mu2.setup(sphereDataConfig, 2);
		sphSolver_inv_one_minus_mu2.solver_component_rexi_z1(1.0, 1.0);	// (1.0
		sphSolver_inv_one_minus_mu2.solver_component_rexi_z2(-1.0, 1.0);	//      - mu^2)


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
#pragma omp parallel for schedule(guided)
#endif
		for (int m = i_sph_data.sphereDataConfig->spectral_modes_m_max; m >= 0; m--)
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
		return inv_one_minus_mu2(spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data));
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
	 * Divergence Operator along longitude for Robert function formulation
	 *
	 * This computes
	 * 	1/cos^2(phi)  d/dlambda U
	 */
	SphereData robert_div_lon(
			const SphereData &i_sph_data
	)
	{
		return inv_one_minus_mu2(diff_lon(i_sph_data));
	}


	SphereData inv_one_minus_mu2(
			const SphereData &i_sph_data
	)
	{
#if 1

		return sphSolver_inv_one_minus_mu2.solve(i_sph_data);

#elif 1
		i_sph_data.request_data_spectral();


		/*
		 * Scale in physical space
		 * => This leads to spurious modes!
		 */

		SphereData out(i_sph_data);

		// Physical space
		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, double &o_data)
				{
					o_data /= cos_phi*cos_phi;
				}
			);

		return out;

#else
		i_sph_data.request_data_spectral();


		// TODO: This is wrong!!!
		SphereData out(i_sph_data.sphereDataConfig);

		// Physical space
		out.spectral_update_lambda(
			[&i_sph_data](
					int n, int m,
					std::complex<double> &o_data
			)
			{
				if (m == 0)
				{
					o_data = 0;
					return;
				}

				o_data = -1.0/(2.0*(double)m);

				double A = std::sqrt(
						(n+m+2.0)*(n+m+1.0)*(2.0*n+1.0)/(2.0*n+3.0)
					);
				double B = std::sqrt(
						(n-m+2.0)*(n-m+1.0)*(2.0*n+1.0)/(2.0*n+3.0)
					);

				o_data *= 	A*i_sph_data.spectral_get(n+1, m+1)
							+ B*i_sph_data.spectral_get(n+1, m-1);
			}
		);

		return out;
#endif
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
		 *   1/(1-sin^2(phi)) * cos^2(phi) * d/d mu f(lambda,mu)
		 */
		return inv_one_minus_mu2(spectral_cosphi2_diff_lat_mu(i_sph_data));
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
		return inv_one_minus_mu2(robert_grad_lon(i_sph_data));
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
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 */
	SphereData robert_grad_lat_M(
			const SphereData &i_sph_data
	)
	{
		return inv_one_minus_mu2(robert_grad_lat(i_sph_data));
	}


	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 */
	SphereData robert_grad_M(
			const SphereData &i_phi,
			const SphereData &i_u,
			const SphereData &i_v
	)
	{
		return inv_one_minus_mu2(
				diff_lon(i_phi)*i_u +
				spectral_one_minus_mu_squared_diff_lat_mu(i_phi)*i_v
			);
	}






	SphereData spec_one_minus_sinphi_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	SphereData spectral_cosphi2_diff_lat_mu(
			const SphereData &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * (1-mu^2) d/dmu ()
	 */
	SphereData spectral_one_minus_mu_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)
	{
		i_sph_data.request_data_spectral();
		const SphereDataConfig *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereData out_sph_data = SphereData(sphereDataConfig);

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
						((-n+1.0)*R(n-1,m))*i_sph_data.spectral_get(n-1, m) +
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
	)	const
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
		const SphereDataConfig *sphereDataConfig = i_sph_data.sphereDataConfig;
		i_sph_data.request_data_spectral();

		SphereData out_sph_data = SphereData(sphereDataConfig);


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
		SphereData out_sph_data = spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);

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
	SphereData div_lat(
			const SphereData &i_sph_data
	)
	{
		SphereData out_sph_data(i_sph_data);



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

		return out_sph_data;
	}


	/**
	 * Laplace operator
	 */
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
#if 1

		return robert_div_lon(i_lon) + robert_div_lat(i_lat);

#else

		return inv_one_minus_mu2(
				diff_lon(i_lon)
				+ spectral_cosphi2_diff_lat_mu(i_lat)
			);
#endif
	}
};



#endif /* SPHOPERATORS_HPP_ */
