/*
 * SPHOperators.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SPHOPERATORS_HPP_
#define SPHOPERATORS_HPP_

#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/sphere/SphereHelpers_SPHIdentities.hpp>
#include <sweet/SimulationVariables.hpp>



class SphereOperators_SphereData	:
	public SphereHelpers_SPHIdentities
{
	friend SphereData_Config;

public:
	const SphereData_Config *sphereDataConfig;

	// Coriolis effect in physical space
	SphereData_Physical fg;


private:
	double r;
	double ir;


public:
	SphereOperators_SphereData(
		SphereData_Config *i_sphereDataConfig,
		const SimulationVariables::SimulationCoefficients *i_simCoeffs
	)
	{
		setup(i_sphereDataConfig, i_simCoeffs);
	}


public:
	SphereOperators_SphereData()
	{
	}



public:
	void setup(
		const SphereData_Config *i_sphereDataConfig,
		const SimulationVariables::SimulationCoefficients *i_simCoeffs
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		r = i_simCoeffs->sphere_radius;
		ir = 1.0/r;

		fg.setup_if_required(i_sphereDataConfig);
		if (i_simCoeffs->sphere_use_fsphere)
		{
			fg.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = i_simCoeffs->sphere_fsphere_f0;
				}
			);
		}
		else
		{
			fg.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = mu*2.0*i_simCoeffs->sphere_rotating_coriolis_omega;
				}
			);
		}

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




	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void robert_vortdiv_to_uv(
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v

	)	const
	{
		assert(o_u.physical_space_data != nullptr);
		assert(o_v.physical_space_data != nullptr);

		SphereData_Spectral psi = inv_laplace(i_vrt)*ir;
		SphereData_Spectral chi = inv_laplace(i_div)*ir;

		o_u.setup_if_required(i_vrt.sphereDataConfig);
		o_v.setup_if_required(i_vrt.sphereDataConfig);
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
			const SphereData_Spectral &i_phi,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v,
			double i_radius

	)	const
	{
		assert(o_u.physical_space_data != nullptr);
		assert(o_v.physical_space_data != nullptr);

		double ir = 1.0/i_radius;

		SphereData_Spectral psi(sphereDataConfig);
		psi.spectral_set_zero();


		o_u.setup_if_required(i_phi.sphereDataConfig);
		o_v.setup_if_required(i_phi.sphereDataConfig);
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
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v

	)	const
	{
		assert(o_u.physical_space_data != nullptr);
		assert(o_v.physical_space_data != nullptr);

		SphereData_Spectral psi = inv_laplace(i_vrt)*ir;
		SphereData_Spectral chi = inv_laplace(i_div)*ir;

		#if SWEET_DEBUG
			#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
				if (omp_in_parallel())
					FatalError("IN PARALLEL REGION!!!");
			#endif
		#endif


		o_u.setup_if_required(i_vrt.sphereDataConfig);
		o_v.setup_if_required(i_vrt.sphereDataConfig);

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



	void vortdiv_to_uv(
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v,
			bool i_robert_formulation

	)	const
	{
		assert(o_u.physical_space_data != nullptr);
		assert(o_v.physical_space_data != nullptr);

		o_u.setup_if_required(i_vrt.sphereDataConfig);
		o_v.setup_if_required(i_vrt.sphereDataConfig);

		if (!i_robert_formulation)
		{
			vortdiv_to_uv(
					i_vrt, i_div,
					o_u, o_v
			);
		}
		else
		{
			robert_vortdiv_to_uv(
					i_vrt, i_div,
					o_u, o_v
			);
		}
	}



	/**
	 * Convert spectral scalar field to physical one
	 */
	void scalar_spectral_to_physical(
			const SphereData_Spectral &i_spectral,
			SphereData_Physical &o_physical

	)	const
	{
		assert(o_physical.physical_space_data != nullptr);

		o_physical.setup_if_required(i_spectral.sphereDataConfig);

		SH_to_spat(
				sphereDataConfig->shtns,
				i_spectral.spectral_space_data,
				o_physical.physical_space_data
		);
	}



	/**
	 * Convert spectral scalar field to physical one
	 */
	SphereData_Physical scalar_spectral_to_physical(
			const SphereData_Spectral &i_spectral

	)	const
	{
		SphereData_Physical retdata(i_spectral.sphereDataConfig);

		SH_to_spat(
				sphereDataConfig->shtns,
				i_spectral.spectral_space_data,
				retdata.physical_space_data
		);

		return retdata;
	}



	/**
	 * Convert spectral scalar field to physical one
	 */
	void scalar_physical_to_spectral(
			const SphereData_Physical &i_physical,
			SphereData_Spectral &o_spectral
	)	const
	{
		assert(o_spectral.spectral_space_data != nullptr);

		o_spectral.setup_if_required(i_physical.sphereDataConfig);

		spat_to_SH(
				sphereDataConfig->shtns,
				i_physical.physical_space_data,
				o_spectral.spectral_space_data
		);
	}



	/**
	 * Convert spectral scalar field to physical one
	 */
	SphereData_Spectral scalar_physical_to_spectral(
			const SphereData_Physical &i_physical
	)	const
	{
		SphereData_Spectral retdata(i_physical.sphereDataConfig);

		spat_to_SH(
				sphereDataConfig->shtns,
				i_physical.physical_space_data,
				retdata.spectral_space_data
		);

		return retdata;
	}



	SphereData_Spectral robert_uv_to_vort(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v

	)	const
	{
		SphereData_Spectral tmp(sphereDataConfig);
		SphereData_Spectral vort(sphereDataConfig);

		SphereData_Physical ug = i_u;
		SphereData_Physical vg = i_v;

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



	SphereData_Spectral uv_to_vort(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v

	)	const
	{
		SphereData_Spectral tmp(sphereDataConfig);
		SphereData_Spectral vort(sphereDataConfig);

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



	SphereData_Spectral uv_to_div(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v
	)	const
	{
		SphereData_Spectral tmp(sphereDataConfig);
		SphereData_Spectral div(sphereDataConfig);

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
				tmp.spectral_space_data,
				div.spectral_space_data
		);
		shtns_robert_form(sphereDataConfig->shtns, 1);

		return laplace(div)*r;
	}



	void robert_uv_to_vortdiv(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div

	)	const
	{
		assert(o_vrt.spectral_space_data != nullptr);
		assert(o_div.spectral_space_data != nullptr);

		// create copy since the data is modified!
		SphereData_Physical ug = i_u;
		SphereData_Physical vg = i_v;

		o_vrt.setup_if_required(i_u.sphereDataConfig);
		o_div.setup_if_required(i_u.sphereDataConfig);

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				ug.physical_space_data,
				vg.physical_space_data,
				o_vrt.spectral_space_data,
				o_div.spectral_space_data
		);

		o_vrt = laplace(o_vrt)*r;
		o_div = laplace(o_div)*r;
	}



	void uv_to_vortdiv(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div

	)	const
	{
		assert(o_vrt.spectral_space_data != nullptr);
		assert(o_div.spectral_space_data != nullptr);


		o_vrt.setup_if_required(i_u.sphereDataConfig);
		o_div.setup_if_required(i_u.sphereDataConfig);

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
				o_vrt.spectral_space_data,
				o_div.spectral_space_data
		);
		shtns_robert_form(sphereDataConfig->shtns, 1);

		o_vrt = laplace(o_vrt)*r;
		o_div = laplace(o_div)*r;
	}



	void uv_to_vortdiv(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div,
			bool i_robert_formulation

	)	const
	{
		assert(o_vrt.spectral_space_data != nullptr);
		assert(o_div.spectral_space_data != nullptr);

		o_vrt.setup_if_required(i_u.sphereDataConfig);
		o_div.setup_if_required(i_u.sphereDataConfig);

		if (!i_robert_formulation)
		{
			uv_to_vortdiv(
					i_u, i_v,
					o_vrt, o_div
			);
		}
		else
		{
			robert_uv_to_vortdiv(
					i_u, i_v,
					o_vrt, o_div
			);
		}
	}



	SphereData_Spectral spectral_one_minus_sinphi_squared_diff_lat_mu(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	SphereData_Spectral spectral_cosphi2_diff_lat_mu(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	/**
	 * (1-mu^2) d/dmu ()
	 */
	SphereData_Spectral spectral_one_minus_mu_squared_diff_lat_mu(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		const SphereData_Config *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereData_Spectral out_sph_data(sphereDataConfig);

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
	SphereData_Spectral mu(
			const SphereData_Spectral &i_sphere_data
	)	const
	{
		const SphereData_Config *sphereDataConfig = i_sphere_data.sphereDataConfig;

		SphereData_Spectral out_sph_data = SphereData_Spectral(sphereDataConfig);


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
	 * mu*mu*F(\lambda,\mu)
	 */
	SphereData_Spectral mu2(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		const SphereData_Config *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereData_Spectral out_sph_data = SphereData_Spectral(sphereDataConfig);


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
	SphereData_Spectral laplace(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		SphereData_Spectral out_sph_data(i_sph_data);

		out_sph_data.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					o_data *= -(double)n*((double)n+1.0)*(ir*ir);
				}
			);

		return out_sph_data;
	}



	/**
	 * Laplace operator
	 */
	SphereData_Spectral inv_laplace(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		SphereData_Spectral out(i_sph_data);

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
