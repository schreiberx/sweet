/*
 * SPHOperators.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <M.Schreiber@exeter.ac.uk>
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


	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolver_inv_one_minus_mu2;

	double r;
	double ir;


public:
	SphereOperators(
		SphereDataConfig *i_sphereDataConfig,
		double i_earth_radius
	)
	{
		setup(i_sphereDataConfig, i_earth_radius);
	}


public:
	SphereOperators()
	{
	}



public:
	void setup(
		SphereDataConfig *i_sphereDataConfig,
		double i_earth_radius
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		sphSolver_inv_one_minus_mu2.setup(sphereDataConfig, 2);
		sphSolver_inv_one_minus_mu2.solver_component_rexi_z1(1.0, 1.0);	// (1.0
		sphSolver_inv_one_minus_mu2.solver_component_rexi_z2(-1.0, 1.0);	//      - mu^2)

		r = i_earth_radius;
		ir = 1.0/r;

#if 1
		double *mx = new double[2*sphereDataConfig->shtns->nlm];
		st_dt_matrix(sphereDataConfig->shtns, mx);

		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				double a = (-n+1.0)*R(n-1,m);
				double b = (n+2.0)*S(n+1,m);

				if (n+1 > sphereDataConfig->spectral_modes_n_max)
					b = 0;

				//std::cout << idx << ": " << a << "\t" << mx[idx*2+0] << std::endl;
				//std::cout << idx << ": " << b << "\t" << mx[idx*2+1] << std::endl;

				double errora = std::abs(a+mx[idx*2+0]);
				double errorb = std::abs(b+mx[idx*2+1]);

				if (errora > 1e-12 || errorb > 1e-12)
				{
					std::cout << idx << ": n=" << n << ", m=" << m << " | "<< errora << "\t" << errorb << std::endl;
					FatalError("SAFETY CHECK NOT SUCCESSFUL");
				}

				idx++;
			}
		}
		delete [] mx;
#endif

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
	)	const
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
	)	const
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
	)	const
	{
		return grad_lat(i_sph_data);
	}



	/**
	 * Compute gradient component along longitude (lambda)
	 */
	SphereData grad(
			const SphereData &i_phi,
			const SphereData &i_u,
			const SphereData &i_v
	)	const
	{
		return i_u*grad_lon(i_phi) + i_v*grad_lat(i_phi);
	}



	/**
	 * Compute gradient component along longitude (lambda)
	 *
	 * 1.0/sqrt(1-mu*mu) d/dlambda()
	 */
	SphereData grad_lon(
			const SphereData &i_sph_data
	)	const
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
	 *
	 * 	= 1/(1-mu^2) d/dlambda U
	 */
	SphereData robert_div_lon(
			const SphereData &i_sph_data
	)	const
	{
		return inv_one_minus_mu2(diff_lon(i_sph_data));
	}



	/**
	 * Multiply with cos(phi)
	 */
	SphereData toRobert(
			const SphereData &i_sph_data
	)	const
	{
		SphereData out(i_sph_data);

		// Physical space
		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, double &o_data)
				{
					o_data *= cos_phi;
				}
			);

		return out;
	}


	/**
	 * Divide by cos(phi)
	 */
	SphereData fromRobert(
			const SphereData &i_sph_data
	)	const
	{
		SphereData out(i_sph_data);

		// Physical space
		out.physical_update_lambda_cosphi_grid(
				[](double lambda, double cos_phi, double &o_data)
				{
					o_data /= cos_phi;
				}
			);

		return out;
	}


	SphereData inv_one_minus_mu2(
			const SphereData &i_sph_data
	)	const
	{
#if 1

		return sphSolver_inv_one_minus_mu2.solve(i_sph_data);

#elif 1
		/*
		 * WARNING: THIS VERSION RESULTS IN REAL EIGENVALUES!!!!
		 */
//		i_sph_data.request_data_spectral();


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
#endif
	}



	void uv_to_stream_potential(
			const SphereData &i_u,
			const SphereData &i_v,
			SphereData &o_stream,
			SphereData &o_potential

	)	const
	{
		i_u.request_data_physical();
		i_v.request_data_physical();

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				o_stream.spectral_space_data,
				o_potential.spectral_space_data
		);

		o_stream.physical_space_data_valid = false;
		o_stream.spectral_space_data_valid = true;

		o_potential.physical_space_data_valid = false;
		o_potential.spectral_space_data_valid = true;
	}



	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void robert_vortdiv_to_uv(
			const SphereData &i_vrt,
			const SphereData &i_div,
			SphereDataPhysical &o_u,
			SphereDataPhysical &o_v

	)	const
	{
		i_vrt.request_data_spectral();
		i_div.request_data_spectral();

		SphereData psi = inv_laplace(i_vrt)*ir;
		SphereData chi = inv_laplace(i_div)*ir;

		SHsphtor_to_spat(
				sphereDataConfig->shtns,
				psi.spectral_space_data,
				chi.spectral_space_data,
				o_u.physical_space_data,
				o_v.physical_space_data
		);

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

	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void vortdiv_to_uv(
			const SphereData &i_vrt,
			const SphereData &i_div,
			SphereDataPhysical &o_u,
			SphereDataPhysical &o_v

	)	const
	{

		i_vrt.request_data_spectral();
		i_div.request_data_spectral();

		SphereData psi = inv_laplace(i_vrt)*ir;
		SphereData chi = inv_laplace(i_div)*ir;

		SHsphtor_to_spat(
				sphereDataConfig->shtns,
				psi.spectral_space_data,
				chi.spectral_space_data,
				o_u.physical_space_data,
				o_v.physical_space_data
		);
	}



	SphereData robert_uv_to_vort(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v

	)	const
	{
		SphereDataPhysical ug = i_u;

		ug.physical_update_lambda_cosphi_grid(
			[&](double lon, double phi, double &o_data)
			{
				o_data /= phi;
			}
		);

		SphereDataPhysical vg = i_v;
		vg.physical_update_lambda_cosphi_grid(
			[&](double lon, double phi, double &o_data)
			{
				o_data /= phi;
			}
		);

		SphereData tmp(sphereDataConfig);
		SphereData vort(sphereDataConfig);

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				ug.physical_space_data,
				vg.physical_space_data,
				vort.spectral_space_data,
				tmp.spectral_space_data
		);

		vort.physical_space_data_valid = false;
		vort.spectral_space_data_valid = true;

		return laplace(vort)*r;
	}



	SphereData uv_to_vort(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v

	)	const
	{
		SphereData tmp(sphereDataConfig);
		SphereData vort(sphereDataConfig);

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				vort.spectral_space_data,
				tmp.spectral_space_data
		);

		vort.physical_space_data_valid = false;
		vort.spectral_space_data_valid = true;

		return laplace(vort)*r;
	}


	void robert_uv_to_vortdiv(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v,
			SphereData &o_vort,
			SphereData &o_div

	)	const
	{
		SphereDataPhysical ug = i_u;

		ug.physical_update_lambda_cosphi_grid(
			[&](double lon, double phi, double &o_data)
			{
				o_data /= phi;
			}
		);

		SphereDataPhysical vg = i_v;
		vg.physical_update_lambda_cosphi_grid(
			[&](double lon, double phi, double &o_data)
			{
				o_data /= phi;
			}
		);

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				ug.physical_space_data,
				vg.physical_space_data,
				o_vort.spectral_space_data,
				o_div.spectral_space_data
		);

		o_vort.physical_space_data_valid = false;
		o_vort.spectral_space_data_valid = true;

		o_div.physical_space_data_valid = false;
		o_div.spectral_space_data_valid = true;

		o_vort = laplace(o_vort)*r;
		o_div = laplace(o_div)*r;
	}




	void uv_to_vortdiv(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v,
			SphereData &o_stream,
			SphereData &o_potential

	)	const
	{
		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				o_stream.spectral_space_data,
				o_potential.spectral_space_data
		);

		o_stream.physical_space_data_valid = false;
		o_stream.spectral_space_data_valid = true;

		o_potential.physical_space_data_valid = false;
		o_potential.spectral_space_data_valid = true;

		o_stream = laplace(o_stream)*r;
		o_potential = laplace(o_potential)*r;
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
	)	const
	{
		return diff_lon(i_sph_data);
	}



	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 *
	 * 1.0/(1-mu*mu) d/dlambda Phi
	 */
	SphereData robert_grad_lon_M(
			const SphereData &i_sph_data
	)	const
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
	)	const
	{
		// Entirely in spectral space
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 *
	 * d/dmu Phi
	 */
	SphereData robert_grad_lat_M(
			const SphereData &i_sph_data
	)	const
	{
		return inv_one_minus_mu2(robert_grad_lat(i_sph_data));
	}


	/**
	 */
	SphereData robert_grad(
			const SphereData &i_phi,
			const SphereData &i_u,
			const SphereData &i_v
	)	const
	{
		return
				diff_lon(i_phi)*i_u +
				spectral_one_minus_mu_squared_diff_lat_mu(i_phi)*i_v
			;
	}

	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 */
	SphereData robert_grad_M(
			const SphereData &i_phi,
			const SphereData &i_u,
			const SphereData &i_v
	)	const
	{
		return inv_one_minus_mu2(
				diff_lon(i_phi)*i_u +
				spectral_one_minus_mu_squared_diff_lat_mu(i_phi)*i_v
			);
	}



	/**
	 * A function following Richies paper
	 */
	SphereData ritchie_A(
			const SphereData &i_u,
			const SphereData &i_v,
			const SphereData &i_phi
	)	const
	{
		return inv_one_minus_mu2(
				(
						i_u*diff_lon(i_phi) +
						i_v*spectral_one_minus_mu_squared_diff_lat_mu(i_phi)
				)
			);
	}






	SphereData spectral_one_minus_sinphi_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)	const
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	SphereData spectral_cosphi2_diff_lat_mu(
			const SphereData &i_sph_data
	)	const
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * (1-mu^2) d/dmu ()
	 */
	SphereData spectral_one_minus_mu_squared_diff_lat_mu(
			const SphereData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();
		const SphereDataConfig *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereData out_sph_data(sphereDataConfig);

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
						((double)(-n+1.0)*R(n-1,m))*i_sph_data.spectral_get(n-1, m) +
						((double)(n+2.0)*S(n+1,m))*i_sph_data.spectral_get(n+1, m);

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
	)	const
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
	)	const
	{
		SphereData out_sph_data(i_sph_data);
		out_sph_data.request_data_spectral();

		out_sph_data.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					o_data *= -(double)n*((double)n+1.0)*ir*ir;
				}
			);

		return out_sph_data;
	}


	/**
	 * Laplace operator
	 */
	SphereData inv_laplace(
			const SphereData &i_sph_data
	)	const
	{
		SphereData out(i_sph_data);

		out.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					if (n != 0)
						o_data /= -(double)n*((double)n+1.0)*ir*ir;
					else
						o_data = 0;
				}
			);

		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

		return out;
	}

#if 1
public:
	/**
	 * Compute vorticity
	 *
	 * \eta = div_lon(V_lat) - div_lat(V_lon)
	 */
	SphereData vort(
			const SphereData &i_lon,
			const SphereData &i_lat
	)	const
	{
		return div_lon(i_lat) - div_lat(i_lon);
	}


public:
	/**
	 * Compute vorticity
	 *
	 * \eta = robert_div_lon(V_lat) - robert_div_lat(V_lon)
	 */
	SphereData robert_vort(
			const SphereData &i_lon,
			const SphereData &i_lat
	)	const
	{
		return robert_div_lon(i_lat) - robert_div_lat(i_lon);
	}
#endif

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


public:
	/**
	 * Compute divergence
	 *
	 * \delta = robert_div_lon(i_lon) + robert_div_lan(i_lan)
	 */
	SphereData robert_div(
			const SphereData &i_lon,
			const SphereData &i_lat
	)	const
	{
#if 0

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
