/*
 * SPHOperators.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SPHOPERATORS_HPP_
#define SPHOPERATORS_HPP_

#include <sweet/MemBlockAlloc.hpp>
#include "../sphere/SphereDataSpectral.hpp"
#include "../sphere/SphereSPHIdentities.hpp"
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>


class SphereOperators	:
	public SphereSPHIdentities
{
	friend SphereDataConfig;

public:
	const SphereDataConfig *sphereDataConfig;

private:

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
		const SphereDataConfig *i_sphereDataConfig,
		double i_earth_radius
	)
	{
		sphereDataConfig = i_sphereDataConfig;

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

	}



#if 0

public:
	/**
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
	SphereDataSpectral diff_lon(
			const SphereDataSpectral &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SphereDataSpectral out_sph_data(i_sph_data.sphereDataConfig);

		// compute d/dlambda in spectral space
		SWEET_THREADING_SPACE_PARALLEL_FOR
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
#endif



	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void robert_vortdiv_to_uv(
			const SphereDataSpectral &i_vrt,
			const SphereDataSpectral &i_div,
			SphereDataPhysical &o_u,
			SphereDataPhysical &o_v

	)	const
	{
		SphereDataSpectral psi = inv_laplace(i_vrt)*ir;
		SphereDataSpectral chi = inv_laplace(i_div)*ir;

		SHsphtor_to_spat(
				sphereDataConfig->shtns,
				psi.spectral_space_data,
				chi.spectral_space_data,
				o_u.physical_space_data,
				o_v.physical_space_data
		);
	}


	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void robert_grad_to_vec(
			const SphereDataSpectral &i_phi,
			SphereDataPhysical &o_u,
			SphereDataPhysical &o_v,
			double i_radius

	)	const
	{
		double ir = 1.0/i_radius;

		SphereDataSpectral psi(sphereDataConfig);
		psi.spectral_set_zero();

		SphereDataPhysical u(sphereDataConfig);
		SphereDataPhysical v(sphereDataConfig);

		SHsphtor_to_spat(
						sphereDataConfig->shtns,
						psi.spectral_space_data,
						i_phi.spectral_space_data,
						o_u.physical_space_data,
						o_v.physical_space_data
				);

		o_u *= ir;
		o_v *= ir;
	}



	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void vortdiv_to_uv(
			const SphereDataSpectral &i_vrt,
			const SphereDataSpectral &i_div,
			SphereDataPhysical &o_u,
			SphereDataPhysical &o_v

	)	const
	{
		SphereDataSpectral psi = inv_laplace(i_vrt)*ir;
		SphereDataSpectral chi = inv_laplace(i_div)*ir;

		#if SWEET_DEBUG
			#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
				if (omp_in_parallel())
					FatalError("IN PARALLEL REGION!!!");
			#endif
		#endif

		shtns_robert_form(sphereDataConfig->shtns, 0);
		SHsphtor_to_spat(
				sphereDataConfig->shtns,
				psi.spectral_space_data,
				chi.spectral_space_data,
				o_u.physical_space_data,
				o_v.physical_space_data
		);
		shtns_robert_form(sphereDataConfig->shtns, 1);
	}



	SphereDataSpectral robert_uv_to_vort(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v

	)	const
	{
		SphereDataSpectral tmp(sphereDataConfig);
		SphereDataSpectral vort(sphereDataConfig);

		SphereDataPhysical ug = i_u;
		SphereDataPhysical vg = i_v;

		shtns_robert_form(sphereDataConfig->shtns, 1);
		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				ug.physical_space_data,
				vg.physical_space_data,
				vort.spectral_space_data,
				tmp.spectral_space_data
		);

		return laplace(vort)*r;
	}



	SphereDataSpectral uv_to_vort(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v

	)	const
	{
		SphereDataSpectral tmp(sphereDataConfig);
		SphereDataSpectral vort(sphereDataConfig);

		#if SWEET_DEBUG
			#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
				if (omp_in_parallel())
					FatalError("IN PARALLEL REGION!!!");
			#endif
		#endif

		shtns_robert_form(sphereDataConfig->shtns, 0);
		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				vort.spectral_space_data,
				tmp.spectral_space_data
		);
		shtns_robert_form(sphereDataConfig->shtns, 1);

		return laplace(vort)*r;
	}



	void robert_uv_to_vortdiv(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v,
			SphereDataSpectral &o_vort,
			SphereDataSpectral &o_div

	)	const
	{
		SphereDataPhysical ug = i_u;
		SphereDataPhysical vg = i_v;

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				ug.physical_space_data,
				vg.physical_space_data,
				o_vort.spectral_space_data,
				o_div.spectral_space_data
		);

		o_vort = laplace(o_vort)*r;
		o_div = laplace(o_div)*r;
	}




	void uv_to_vortdiv(
			const SphereDataPhysical &i_u,
			const SphereDataPhysical &i_v,
			SphereDataSpectral &o_stream,
			SphereDataSpectral &o_potential

	)	const
	{

		#if SWEET_DEBUG
			#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
				if (omp_in_parallel())
					FatalError("IN PARALLEL REGION!!!");
			#endif
		#endif

		shtns_robert_form(sphereDataConfig->shtns, 0);
		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				o_stream.spectral_space_data,
				o_potential.spectral_space_data
		);
		shtns_robert_form(sphereDataConfig->shtns, 1);

		o_stream = laplace(o_stream)*r;
		o_potential = laplace(o_potential)*r;
	}




	SphereDataSpectral spectral_one_minus_sinphi_squared_diff_lat_mu(
			const SphereDataSpectral &i_sph_data
	)	const
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	SphereDataSpectral spectral_cosphi2_diff_lat_mu(
			const SphereDataSpectral &i_sph_data
	)	const
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * (1-mu^2) d/dmu ()
	 */
	SphereDataSpectral spectral_one_minus_mu_squared_diff_lat_mu(
			const SphereDataSpectral &i_sph_data
	)	const
	{
		const SphereDataConfig *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereDataSpectral out_sph_data(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sph_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] =
						((double)(-n+1.0)*R(n-1,m))*i_sph_data.spectral_get(n-1, m) +
						((double)(n+2.0)*S(n+1,m))*i_sph_data.spectral_get(n+1, m);

				idx++;
			}
		}

		return out_sph_data;
	}



	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SphereDataSpectral mu(
			const SphereDataSpectral &i_sphere_data
	)	const
	{
		const SphereDataConfig *sphereDataConfig = i_sphere_data.sphereDataConfig;

		SphereDataSpectral out_sph_data = SphereDataSpectral(sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] =
							R(n-1,m)*i_sphere_data.spectral_get(n-1, m)
							+ S(n+1,m)*i_sphere_data.spectral_get(n+1, m);

				idx++;
			}
		}

		return out_sph_data;
	}

	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SphereDataSpectral mu2(
			const SphereDataSpectral &i_sph_data
	)	const
	{
		const SphereDataConfig *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereDataSpectral out_sph_data = SphereDataSpectral(sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
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

		return out_sph_data;
	}



	/**
	 * Laplace operator
	 */
	SphereDataSpectral laplace(
			const SphereDataSpectral &i_sph_data
	)	const
	{
		SphereDataSpectral out_sph_data(i_sph_data);

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
	SphereDataSpectral inv_laplace(
			const SphereDataSpectral &i_sph_data
	)	const
	{
		SphereDataSpectral out(i_sph_data);

		out.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					if (n != 0)
						o_data /= -(double)n*((double)n+1.0)*ir*ir;
					else
						o_data = 0;
				}
			);

		return out;
	}


};



#endif /* SPHOPERATORS_HPP_ */
