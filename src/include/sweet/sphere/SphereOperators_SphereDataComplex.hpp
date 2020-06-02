/*
 * SPHOperatorsComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SPHEREOPERATORS_COMPLEX_HPP_
#define SPHEREOPERATORS_COMPLEX_HPP_

#include <sweet/sphere/Convert_SphereDataPhysicalComplex_to_SphereDataPhysical.hpp>
#include <libmath/shtns_inc.hpp>
#include <sweet/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/sphere/SphereHelpers_SPHIdentities.hpp>
#include <sweet/SimulationVariables.hpp>



class SphereOperators_SphereDataComplex	:
		public SphereHelpers_SPHIdentities
{
	friend SphereData_Config;

	const SphereData_Config *sphereDataConfig;

private:
	double r;
	double ir;


	/**
	 * Constructor
	 */
public:
	SphereOperators_SphereDataComplex(
			SphereData_Config *i_sphereDataConfig,
			const SimulationVariables::SimulationCoefficients *i_simCoeffs
	)
	{
		setup(i_sphereDataConfig, i_simCoeffs);
	}


public:
	SphereOperators_SphereDataComplex()	:
		sphereDataConfig(nullptr)
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
	}


public:


	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (I + b D^2) x = rhs
	 */
	inline
	SphereData_SpectralComplex implicit_helmholtz(
			const SphereData_SpectralComplex &i_sphere_data,
			const std::complex<double> &i_b,
			double i_radius
	)	const
	{
		SphereData_SpectralComplex out(i_sphere_data);

		const std::complex<double> b = i_b/(i_radius*i_radius);

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				io_data /= (1.0 + (-b*(double)n*((double)n+1.0)));
			}
		);

		return out;
	}


	/**
	 * Compute multiplication with "J" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_SpectralComplex implicit_J(
			const SphereData_SpectralComplex &i_sphere_data,
			const std::complex<double>& i_dt_two_omega
	)	const
	{
		SphereData_SpectralComplex out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
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
	SphereData_SpectralComplex implicit_Jinv(
			const SphereData_SpectralComplex &i_sphere_data,
			const std::complex<double>& i_dt_two_omega
	)	const
	{
		SphereData_SpectralComplex out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
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
	SphereData_SpectralComplex implicit_F(
			const SphereData_SpectralComplex &i_sphere_data,
			const std::complex<double>& i_dt_two_omega
	)	const
	{
		SphereData_SpectralComplex out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				//std::cout << "* START *******************************************" << std::endl;
				//std::cout << "implicit_F " << n << ", " << m << std::endl;

				out_sph_data.spectral_space_data[idx] = 0;

				// out of boundary check for P(n-1, m)
				if (n-1 >= std::abs(m))
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_minus(n, m)
							* i_sphere_data.spectral_get_(n-1, m);
					//std::cout << "a: x" << std::endl;
				}

				// out of boundary check for P(n+1, m)
				if (n+1 <= i_sphere_data.sphereDataConfig->spectral_modes_n_max)
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_plus(n, m)
							* i_sphere_data.spectral_get_(n+1, m);

					//std::cout << "b: x" << std::endl;
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
	SphereData_SpectralComplex implicit_FJinv(
			const SphereData_SpectralComplex &i_sphere_data,
			const std::complex<double>& i_dt_two_omega
	)	const
	{
		SphereData_SpectralComplex out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				//std::cout << "* START *******************************************" << std::endl;
				//std::cout << "implicit_FJinv " << n << ", " << m << std::endl;

				out_sph_data.spectral_space_data[idx] = 0;

				// Out of boundary check for P(n-1, m)
				if (n-1 >= std::abs(m))
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_minus(n, m)
							/ implicit_J_scalar(n-1, m, i_dt_two_omega)
							* i_sphere_data.spectral_get_(n-1, m);

					//std::cout << "a: x" << std::endl;
				}

				// Out of boundary check for P(n+1, m)
				if (n+1 <= i_sphere_data.sphereDataConfig->spectral_modes_n_max)
				{
					out_sph_data.spectral_space_data[idx] +=
							i_dt_two_omega
							* implicit_f_plus(n, m)
							/ implicit_J_scalar(n+1, m, i_dt_two_omega)
							* i_sphere_data.spectral_get_(n+1, m);

					//std::cout << "b: x" << std::endl;
				}

				idx++;
			}
		}

		return out_sph_data;
	}



	/**
	 * Compute multiplication with "L" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_SpectralComplex implicit_L(
			const SphereData_SpectralComplex &i_sphere_data,
			const std::complex<double>& i_dt
	)	const
	{
		SphereData_SpectralComplex out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out_sph_data[idx] = i_dt*ir*ir*Tcomplex(n*(n+1)) * i_sphere_data[idx];
				idx++;
			}
		}

		return out_sph_data;
	}



	/**
	 * Compute multiplication with "Linv" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
	SphereData_SpectralComplex implicit_Linv(
			const SphereData_SpectralComplex &i_sphere_data,
			const std::complex<double>& i_dt
	)	const
	{
		SphereData_SpectralComplex out_sph_data(i_sphere_data.sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sphere_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sphere_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				if (n == 0)
					out_sph_data[idx] = 0.0;
				else
					out_sph_data[idx] = 1.0/(i_dt*ir*ir*Treal(n*(n+1))) * i_sphere_data[idx];

				idx++;
			}
		}

		return out_sph_data;
	}

	/**
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
	SphereData_SpectralComplex diff_lon(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		SphereData_SpectralComplex out(i_sph_data.sphereDataConfig);

		// compute d/dlambda in spectral space
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[idx] = i_sph_data.spectral_space_data[idx]*std::complex<double>(0, m);
				idx++;
			}
		}

		return out;
	}


	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SphereData_SpectralComplex mu(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		const SphereData_Config *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereData_SpectralComplex out = SphereData_SpectralComplex(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[idx] = 0;

				if (n-1 >= std::abs(m))
					out.spectral_space_data[idx] += R(n-1,m)*i_sph_data.spectral_get_(n-1, m);

				if (n+1 <= i_sph_data.sphereDataConfig->spectral_modes_n_max)
					out.spectral_space_data[idx] += S(n+1,m)*i_sph_data.spectral_get_(n+1, m);

				idx++;
			}
		}

		return out;
	}

	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SphereData_SpectralComplex mu2(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		const SphereData_Config *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereData_SpectralComplex out = SphereData_SpectralComplex(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[idx] = 0;

				if (n-2 >= std::abs(m))
					out.spectral_space_data[idx] += A(n-2,m)*i_sph_data.spectral_get_(n-2, m);

				out.spectral_space_data[idx] += B(n+0,m)*i_sph_data.spectral_get_(n+0, m);

				if (n+2 <= i_sph_data.sphereDataConfig->spectral_modes_n_max)
					out.spectral_space_data[idx] += C(n+2,m)*i_sph_data.spectral_get_(n+2, m);

				idx++;
			}
		}

		return out;
	}


#if 0
	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d mu f(lambda,mu)
	 *
	 * sqrt(1-mu*mu)*d/dmu P_n^m = ...
	 */
	SphereData_SpectralComplex diff_lat_mu(
			const SphereData_SpectralComplex &i_sph_data
	)
	{

		return inv_one_minus_mu2(spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data));
	}



	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d phi f(lambda,mu)
	 */
	SphereData_SpectralComplex diff_lat_phi(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		return grad_lat(i_sph_data);
	}
#endif

#if 0
	/**
	 * Compute gradient component along longitude (lambda)
	 */
	SphereData_SpectralComplex grad_lon(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		SphereData_SpectralComplex out = diff_lon(i_sph_data);

		out.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					double cos_phi = std::sqrt(1.0-mu*mu);
					o_data /= cos_phi;
				}
		);

		return out;
	}



	SphereData_SpectralComplex inv_one_minus_mu2(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
#if 1

		return sphSolver_inv_one_minus_mu2.solve(i_sph_data);

#else
		SphereData_SpectralComplex out(i_sph_data);

		out.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					o_data /= (1.0-mu*mu);
				}
		);
		return out;
#endif
	}


	SphereData_SpectralComplex spectral_one_minus_mu_squared_diff_lat_mu(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		i_sph_data.request_data_spectral();
		const SphereData_Config *sphereDataConfig = i_sph_data.sphereDataConfig;

		SphereData_SpectralComplex out(sphereDataConfig);


		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				/**
				 * TODO: Optimize me!
				 */
				out.spectral_space_data[idx] =
						((-n+1.0)*R(n-1,m))*i_sph_data.spectral_get_DEPRECATED(n-1, m) +
						((n+2.0)*S(n+1,m))*i_sph_data.spectral_get_DEPRECATED(n+1, m);

				idx++;
			}
		}

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

		return out;
	}


	/**
	 * Compute gradient component along latitude
	 */
	SphereData_SpectralComplex grad_lat(
			const SphereData_SpectralComplex &i_sph_data
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
		SphereData_SpectralComplex out = spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out.physical_update_lambda_gaussian_grid(
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
		out.request_data_physical();
		out.physical_update_lambda_gaussian_grid(
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

		return out;
	}



	/**
	 * Divergence Operator along longitude
	 *
	 * Identical to gradient operator along longitude
	 */
	SphereData_SpectralComplex div_lon(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		return grad_lon(i_sph_data);
	}




	/**
	 * Divergence Operator along latitude
	 *
	 * d(sqrt(1-mu*mu)*F)/dmu
	 */
	SphereData_SpectralComplex div_lat(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		SphereData_SpectralComplex out(i_sph_data);

#if 1
		out.physical_update_lambda_cogaussian_grid(
				[](double lambda, double comu, std::complex<double> &o_data)
				{
					//o_data *= cos(phi);
					o_data *= comu;
				}
			);
#else
		out.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					//o_data *= cos(phi);
					o_data *= std::sqrt(1.0-mu*mu);
				}
			);
#endif
		// grad_lat = diff_lat_phi

#if 1

		out = spectral_one_minus_mu_squared_diff_lat_mu(out);

		out.physical_update_lambda_gaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					//o_data /= mu;
					o_data /= (1.0-mu*mu);
				}
			);

#else
		out = grad_lat(out);

		// undo the sin(theta) which is cos(phi)
		out.physical_update_lambda(
				[](double lambda, double phi, std::complex<double> &o_data)
				{
					//o_data /= mu;
					o_data /= cos(phi);
				}
			);
#endif

		return out;
	}


	inline
	SphereData_SpectralComplex spectral_one_minus_sinphi_squared_diff_lat_mu(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}



	inline
	SphereData_SpectralComplex spectral_cosphi2_diff_lat_mu(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}
#endif


	/**
	 * Laplace operator
	 */
	SphereData_SpectralComplex laplace(
			const SphereData_SpectralComplex &i_sph_data,
			double r = -1
	)	const
	{
		if (r == -1)
			r = this->r;

		double ir = 1.0/r;

		SphereData_SpectralComplex out(i_sph_data);

		out.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					o_data *= -(double)n*((double)n+1.0)*ir*ir;
				}
			);

		return out;
	}


	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void vrtdiv_to_uv(
			const SphereData_SpectralComplex &i_vrt,
			const SphereData_SpectralComplex &i_div,
			SphereData_PhysicalComplex &o_u,
			SphereData_PhysicalComplex &o_v

	)	const
	{
		SphereData_SpectralComplex psi = inv_laplace(i_vrt)*ir;
		SphereData_SpectralComplex chi = inv_laplace(i_div)*ir;

		#if SWEET_DEBUG
			#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
				if (omp_in_parallel())
					SWEETError("IN PARALLEL REGION!!!");
			#endif
		#endif


		o_u.setup_if_required(i_vrt.sphereDataConfig);
		o_v.setup_if_required(i_vrt.sphereDataConfig);

		SHsphtor_to_spat_cplx(
				sphereDataConfig->shtns,
				psi.spectral_space_data,
				chi.spectral_space_data,
				o_u.physical_space_data,
				o_v.physical_space_data
		);
	}


	void uv_to_vrtdiv(
			const SphereData_PhysicalComplex &i_u,
			const SphereData_PhysicalComplex &i_v,
			SphereData_SpectralComplex &o_vrt,
			SphereData_SpectralComplex &o_div

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

		spat_cplx_to_SHsphtor(
				sphereDataConfig->shtns,
				i_u.physical_space_data,
				i_v.physical_space_data,
				o_vrt.spectral_space_data,
				o_div.spectral_space_data
		);

		o_vrt = laplace(o_vrt)*r;
		o_div = laplace(o_div)*r;
	}


#if 0
public:
	/**
	 * Compute vorticity
	 *
	 * \eta = div_lat(V_lon) - div_lon(V_lat)
	 */
	SphereData_SpectralComplex vort(
			const SphereData_SpectralComplex &i_lon,
			const SphereData_SpectralComplex &i_lat
	)
	{
		return div_lon(i_lat) - div_lat(i_lon);
	}


#endif

	/**
	 * Laplace operator
	 */
	SphereData_SpectralComplex inv_laplace(
			const SphereData_SpectralComplex &i_sph_data,
			double i_radius = -1
	)	const
	{
		if (i_radius == -1)
			i_radius = this->r;

		double ir = 1.0/i_radius;

		SphereData_SpectralComplex out(i_sph_data);

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


#if 0

	void uv_to_vrtdiv(
			const SphereData_PhysicalComplex &i_u,
			const SphereData_PhysicalComplex &i_v,
			SphereData_SpectralComplex &o_vort,
			SphereData_SpectralComplex &o_div,

			double r = -1		// FIXME / TODO: Remove me!!!
	)	const
	{
		if (r == -1)
			r = this->r;

		/*
		 * Generate a copy because of destructive SHT operations
		 */
		SphereData_PhysicalComplex ug = i_u;
		SphereData_PhysicalComplex vg = i_v;

		spat_cplx_to_SHsphtor(
				sphereDataConfig->shtns,
				ug.physical_space_data,
				vg.physical_space_data,
				o_vort.spectral_space_data,
				o_div.spectral_space_data
		);

		o_vort = laplace(o_vort, r)*r;
		o_div = laplace(o_div, r)*r;
	}



public:
	/**
	 * Compute divergence
	 *
	 * \delta = div_lon(i_lon) + div_lan(i_lan)
	 */
	SphereData_SpectralComplex div(
			const SphereData_SpectralComplex &i_lon,
			const SphereData_SpectralComplex &i_lat
	)
	{
		return div_lon(i_lon) + div_lat(i_lat);
	}
#endif

};






#endif /* SPHOPERATORS_HPP_ */
