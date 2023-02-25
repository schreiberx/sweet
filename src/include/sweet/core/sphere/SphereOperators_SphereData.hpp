/*
 * SPHOperators.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SPHOPERATORS_HPP_
#define SPHOPERATORS_HPP_

#include <sweet/core/MemBlockAlloc.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereHelpers_SPHIdentities.hpp>
#include <sweet/core/SimulationVariables.hpp>



class SphereOperators_SphereData	:
	public SphereHelpers_SPHIdentities
{
	friend SphereData_Config;
	typedef std::complex<double> Tcomplex;
	typedef double Treal;

public:
	const SphereData_Config *sphereDataConfig;

	// Coriolis effect in physical space
	SphereData_Physical fg;

	// Solely the rotational effect without anything else
	SphereData_Physical mug;


private:
	double r;		// radius
	double r2;		// radius^2

	double ir;		// 1/radius
	double ir2;		// 1/radius^2


public:
	SphereOperators_SphereData(
		SphereData_Config *i_sphereDataConfig,
		const SimulationCoefficients *i_simCoeffs
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
		const SimulationCoefficients *i_simCoeffs
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		r = i_simCoeffs->sphere_radius;
		r2 = r*r;

		ir = 1.0/r;
		ir2 = ir*ir;

		fg.free();
		fg.setup(i_sphereDataConfig);
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

		mug.free();
		mug.setup(i_sphereDataConfig);
		if (i_simCoeffs->sphere_use_fsphere)
		{
			mug.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = i_simCoeffs->sphere_fsphere_f0;
				}
			);
		}
		else
		{
			mug.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = mu;
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
					SWEETError("SAFETY CHECK NOT SUCCESSFUL");
				}

				idx++;
			}
		}
		delete [] mx;
#endif

	}





	/**
	 * Compute gradient
	 */
	void grad_to_vec(
			const SphereData_Spectral &i_phi,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v,
			double i_radius

	)	const
	{
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
	void vrtdiv_to_uv(
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,
			SphereData_Physical &o_u,
			SphereData_Physical &o_v

	)	const
	{
		/* Calculate stream function and velocity potential (multiplied by earth radius) */
		SphereData_Spectral psi = inv_laplace(i_vrt)*ir;
		SphereData_Spectral chi = inv_laplace(i_div)*ir;

		#if SWEET_DEBUG
			#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
				if (omp_in_parallel())
					SWEETError("IN PARALLEL REGION!!!");
			#endif
		#endif


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
	 * Convert spectral scalar field to physical one
	 */
	void scalar_spectral_to_physical(
			const SphereData_Spectral &i_spectral,
			SphereData_Physical &o_physical

	)	const
	{
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
					SWEETError("IN PARALLEL REGION!!!");
			#endif
		#endif

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				vort.spectral_space_data,
				tmp.spectral_space_data
		);

		return laplace(vort)*r;
	}



	/*
	 * Compute Nonlinear advection terms
	 *
	 * U \cdot \grad phi = \div \cdot (V*phi) - \nabla
	 */
	SphereData_Spectral V_dot_grad_scalar(
			const SphereData_Physical &i_u_phys,		///< u velocity
			const SphereData_Physical &i_v_phys,		///< v velocity
			const SphereData_Physical &i_div_phys,		///< divergence in physical space to avoid transformation
			const SphereData_Physical &i_scalar_phys	///< scalar
	)	const
	{
		return uv_to_div(
				i_u_phys*i_scalar_phys,
				i_v_phys*i_scalar_phys
			)
				- i_div_phys*i_scalar_phys;
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
					SWEETError("IN PARALLEL REGION!!!");
			#endif
		#endif

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				tmp.spectral_space_data,
				div.spectral_space_data
		);

		return laplace(div)*r;
	}



	void uv_to_vrtdiv(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div

	)	const
	{
		o_vrt.setup_if_required(i_u.sphereDataConfig);
		o_div.setup_if_required(i_u.sphereDataConfig);

		#if SWEET_DEBUG
			#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
				if (omp_in_parallel())
					SWEETError("IN PARALLEL REGION!!!");
			#endif
		#endif

		spat_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				o_vrt.spectral_space_data,
				o_div.spectral_space_data
		);

		o_vrt = laplace(o_vrt)*r;
		o_div = laplace(o_div)*r;
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
						((double)(-n+1.0)*R(n-1,m))*i_sph_data.spectral_get_DEPRECATED(n-1, m) +
						((double)(n+2.0)*S(n+1,m))*i_sph_data.spectral_get_DEPRECATED(n+1, m);

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
		SphereData_Spectral out_sph_data(i_sphere_data.sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] = 0;

				if (n-1 >= m)
					out_sph_data.spectral_space_data[idx] +=
							R_real(n-1,m)*i_sphere_data.spectral_get_(n-1, m);

				if (n+1 <= i_sphere_data.sphereDataConfig->spectral_modes_n_max)
					out_sph_data.spectral_space_data[idx] +=
							S_real(n+1,m)*i_sphere_data.spectral_get_(n+1, m);

				idx++;
			}
		}

		return out_sph_data;
	}



	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (1 + b D^2) x = rhs
	 */
	inline
	SphereData_Spectral implicit_helmholtz(
			const SphereData_Spectral &i_sphere_data,
			const double &i_b,
			double i_sphere_radius
	)
	{
		SphereData_Spectral out(i_sphere_data);

		const double b = i_b/(i_sphere_radius*i_sphere_radius);

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				/*
				 * Note, the Laplace operator in SPH space is given by
				 * 	-(double)n*((double)n+1.0))/r^2
				 */
				io_data /= (1 + (-b*(double)n*((double)n+1.0)));
			}
		);

		return out;
	}




	/**
	 * Compute multiplication with "J" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_Spectral implicit_J(
			const SphereData_Spectral &i_sphere_data,
			double i_dt_two_omega
	)	const
	{
		SphereData_Spectral out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data[idx] = i_sphere_data[idx] * implicit_J_scalar(n, m, i_dt_two_omega);
				idx++;
			}
		}

		return out_sph_data;
	}


	/**
	 * Compute multiplication with "J^-1" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_Spectral implicit_Jinv(
			const SphereData_Spectral &i_sphere_data,
			double i_dt_two_omega
	)	const
	{
		SphereData_Spectral out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data[idx] =	i_sphere_data.spectral_space_data[idx] / implicit_J_scalar(n, m, i_dt_two_omega);
				idx++;
			}
		}

		return out_sph_data;
	}


	/**
	 * Compute multiplication with "F" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_Spectral implicit_F(
			const SphereData_Spectral &i_sphere_data,
			double i_dt_two_omega
	)	const
	{
		SphereData_Spectral out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] = 0;

				// out of boundary check for P(n-1, m)
				if (n-1 >= m)
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_minus(n, m)
							* i_sphere_data.spectral_get_(n-1, m);
				}

				// out of boundary check for P(n+1, m)
				if (n+1 <= i_sphere_data.sphereDataConfig->spectral_modes_n_max)
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_plus(n, m)
							* i_sphere_data.spectral_get_(n+1, m);
				}

				idx++;
			}
		}

		return out_sph_data;
	}



	/**
	 * Compute multiplication with "F" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_Spectral implicit_FJinv(
			const SphereData_Spectral &i_sphere_data,
			double i_dt_two_omega
	)	const
	{
		SphereData_Spectral out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				out_sph_data.spectral_space_data[idx] = 0;

				// Out of boundary check for P(n-1, m)
				if (n-1 >= m)
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_minus(n, m)
							/ implicit_J_scalar(n-1, m, i_dt_two_omega)
							* i_sphere_data.spectral_get_(n-1, m);
				}

				// Out of boundary check for P(n+1, m)
				if (n+1 <= i_sphere_data.sphereDataConfig->spectral_modes_n_max)
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_plus(n, m)
							/ implicit_J_scalar(n+1, m, i_dt_two_omega)
							* i_sphere_data.spectral_get_(n+1, m);
				}

				idx++;
			}
		}

		return out_sph_data;
	}



	/**
	 * Compute multiplication with "L" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 *
	 * return -dt * D^2 * sphere_data
	 */
	SphereData_Spectral implicit_L(
			const SphereData_Spectral &i_sphere_data,
			double i_dt
	)	const
	{
		SphereData_Spectral out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				/*
				 * Note, the Laplace operator in SPH space is given by
				 * 	-(double)n*((double)n+1.0))/r^2
				 */
				out_sph_data[idx] = i_dt*ir2*(n*(n+1)) * i_sphere_data[idx];
				idx++;
			}
		}

		return out_sph_data;
	}



	/**
	 * Compute multiplication with "Linv" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_Spectral implicit_Linv(
			const SphereData_Spectral &i_sphere_data,
			double i_dt
	)	const
	{
		SphereData_Spectral out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= i_sphere_data.sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes(m, m);

			for (int n = m; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
			{
				if (n == 0)
					out_sph_data[idx] = 0.0;
				else
					out_sph_data[idx] = 1.0/(i_dt*ir2*(n*(n+1))) * i_sphere_data[idx];

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
						+A(n-2,m)*i_sph_data.spectral_get_DEPRECATED(n-2, m)
						+B(n+0,m)*i_sph_data.spectral_get_DEPRECATED(n+0, m)
						+C(n+2,m)*i_sph_data.spectral_get_DEPRECATED(n+2, m)
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
	 * root Laplace operator
	 */
	SphereData_Spectral root_laplace(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		SphereData_Spectral out_sph_data(i_sph_data);

		out_sph_data.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					o_data *= std::sqrt((double)n*((double)n+1.0))*(ir);
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


	/**
	 * inverse root Laplace operator
	 */
	SphereData_Spectral inv_root_laplace(
			const SphereData_Spectral &i_sph_data
	)	const
	{
		SphereData_Spectral out(i_sph_data);

		out.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					if (n != 0)
						o_data /= std::sqrt((double)n*((double)n+1.0))*ir;
					else
						o_data = 0;
				}
			);

		return out;
	}


	/**
	 * Calculates implicit diffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	inline SphereData_Spectral implicit_diffusion(
			const SphereData_Spectral &i_data,
			double i_coef,
			/*int i_order*/
			double i_r
	)
	{
		SphereData_Spectral out = i_data;

		const double scalar = i_coef;
		const double r      = i_r;

		out  = out.spectral_solve_helmholtz(1.0,  -scalar, r);

		return out;
	}

	/**
	 * Calculates implicit hyperdiffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	inline SphereData_Spectral implicit_hyperdiffusion(
			const SphereData_Spectral &i_data,
			double i_coef,
			int i_order,
			double i_r
	)
	{
		SphereData_Spectral out = i_data;

		const double r      = i_r;

		std::array<double, 4> visc_factors;
		if (i_order == 2)
			visc_factors = {-i_coef, 0, 0, 0};
		else if (i_order == 4)
			visc_factors = {0, -i_coef, 0, 0};
		else if (i_order == 6)
			visc_factors = {0, 0, -i_coef, 0};
		else if (i_order == 8)
			visc_factors = {0, 0, 0, -i_coef};
		else
			SWEETError("This viscosity order is not supported: " + std::to_string(i_order));

		out  = out.spectral_solve_helmholtz_higher_order(1.0, visc_factors, r);

		return out;
	}

};



#endif /* SPHOPERATORS_HPP_ */
