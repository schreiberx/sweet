/*
 * SPHOperatorsComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SPHEREOPERATORS_COMPLEX_HPP_
#define SPHEREOPERATORS_COMPLEX_HPP_

#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalComplex.hpp>
#include <sweet/sphere/Convert_SphereDataPhysicalComplex_to_SphereDataPhysical.hpp>
#include <libmath/shtns_inc.hpp>
#include <sweet/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/sphere/SphereHelpers_SPHIdentities.hpp>



class SphereOperators_SphereDataComplex	:
		public SphereHelpers_SPHIdentities
{
	friend SphereData_Config;

	const SphereData_Config *sphereDataConfig;

public:
	SphBandedMatrixPhysicalComplex< std::complex<double> > sphSolver_inv_one_minus_mu2;

private:
	double r;
	double ir;


	/**
	 * Constructor
	 */
public:
	SphereOperators_SphereDataComplex(
			SphereData_Config *i_sphereDataConfig,
			double i_earth_radius
	)
	{
		setup(i_sphereDataConfig, i_earth_radius);
	}


public:
	SphereOperators_SphereDataComplex()	:
		sphereDataConfig(nullptr)
	{

	}


public:
	void setup(
			const SphereData_Config *i_sphereDataConfig,
			double i_earth_radius
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		r = i_earth_radius;
		ir = 1.0/r;

		sphSolver_inv_one_minus_mu2.setup(sphereDataConfig, 2);
		sphSolver_inv_one_minus_mu2.solver_component_rexi_z1(1.0, 1.0);	// (1.0
		sphSolver_inv_one_minus_mu2.solver_component_rexi_z2(-1.0, 1.0);	//      - mu^2)
	}


public:
	/**
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
	SphereData_SpectralComplex diff_lon(
			const SphereData_SpectralComplex &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

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
		out.spectral_space_data_valid = true;
		out.physical_space_data_valid = false;

		return out;
	}



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
						((-n+1.0)*R(n-1,m))*i_sph_data.spectral_get(n-1, m) +
						((n+2.0)*S(n+1,m))*i_sph_data.spectral_get(n+1, m);

				idx++;
			}
		}

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

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
		i_sph_data.request_data_spectral();

		SphereData_SpectralComplex out = SphereData_SpectralComplex(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[idx] =
					R(n-1,m)*i_sph_data.spectral_get(n-1, m)
					+ S(n+1,m)*i_sph_data.spectral_get(n+1, m);

				idx++;
			}
		}

		out.physical_space_data_valid = false;
		out.spectral_space_data_valid = true;

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
		i_sph_data.request_data_spectral();

		SphereData_SpectralComplex out = SphereData_SpectralComplex(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= i_sph_data.sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = i_sph_data.sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[idx] =
						+A(n-2,m)*i_sph_data.spectral_get(n-2, m)
						+B(n+0,m)*i_sph_data.spectral_get(n+0, m)
						+C(n+2,m)*i_sph_data.spectral_get(n+2, m)
						;
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



	/**
	 * Divergence Operator along longitude for robert function formlation
	 *
	 * This computes
	 * 	1/cos^2(phi)  d/dlambda U
	 */
	SphereData_SpectralComplex robert_div_lon(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		return inv_one_minus_mu2(diff_lon(i_sph_data));
	}



	/**
	 * Compute divergence along latitude for robert function formulation
	 *
	 * This computes
	 * 		d/dmu V
	 *
	 * There's no other metric term involved!
	 */
	SphereData_SpectralComplex robert_div_lat(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		/*
		 * Compute
		 *   cos^2(phi) * d/d mu  f(lambda,mu)
		 */
		return inv_one_minus_mu2(spectral_cosphi2_diff_lat_mu(i_sph_data));
	}



	SphereData_SpectralComplex robert_cos2phi_div_lat(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		/*
		 * Compute
		 *   d/d mu  f(lambda,mu)
		 */
		return spectral_cosphi2_diff_lat_mu(i_sph_data);
	}



	/**
	 * Compute gradient component along longitude (lambda) for Robert function formulation
	 *
	 * This computes
	 * 		d/dlambda Phi
	 */
	SphereData_SpectralComplex robert_grad_lon(
			const SphereData_SpectralComplex &i_sph_data
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
	SphereData_SpectralComplex robert_grad_lat(
			const SphereData_SpectralComplex &i_sph_data
	)
	{
		return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
	}


	/**
	 * Special formulation for Robert gradient,
	 * see REXI with spherical harmonics
	 */
	SphereData_SpectralComplex robert_grad_M(
			const SphereData_SpectralComplex &i_phi,
			const SphereData_SpectralComplex &i_u,
			const SphereData_SpectralComplex &i_v
	)
	{
		return inv_one_minus_mu2(
				diff_lon(i_phi)*i_u +
				spectral_one_minus_mu_squared_diff_lat_mu(i_phi)*i_v
			);
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

		i_sph_data.request_data_spectral();

		SphereData_SpectralComplex out(i_sph_data);

		out.spectral_update_lambda(
				[&](int n, int m, std::complex<double> &o_data)
				{
					o_data *= -(double)n*((double)n+1.0)*ir*ir;
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
	SphereData_SpectralComplex vort(
			const SphereData_SpectralComplex &i_lon,
			const SphereData_SpectralComplex &i_lat
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
	SphereData_SpectralComplex robert_vort(
			const SphereData_SpectralComplex &i_lon,
			const SphereData_SpectralComplex &i_lat
	)
	{
		return robert_div_lon(i_lat) - robert_div_lat(i_lon);
	}


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
		out.request_data_spectral();

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



	void robert_uv_to_vortdiv(
			const SphereData_PhysicalComplex &i_u,
			const SphereData_PhysicalComplex &i_v,
			SphereData_SpectralComplex &o_vort,
			SphereData_SpectralComplex &o_div,

			double r = -1
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
		o_vort.spectral_space_data_valid = true;
		o_vort.physical_space_data_valid = false;
		o_div.spectral_space_data_valid = true;
		o_div.physical_space_data_valid = false;

		o_vort = laplace(o_vort)*r;
		o_div = laplace(o_div, r)*r;
	}



	void uv_to_vortdiv(
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

		shtns_robert_form(sphereDataConfig->shtns, 0);
		spat_cplx_to_SHsphtor(
				sphereDataConfig->shtns,
				ug.physical_space_data,
				vg.physical_space_data,
				o_vort.spectral_space_data,
				o_div.spectral_space_data
		);
		shtns_robert_form(sphereDataConfig->shtns, 1);
		o_vort.spectral_space_data_valid = true;
		o_vort.physical_space_data_valid = false;
		o_div.spectral_space_data_valid = true;
		o_div.physical_space_data_valid = false;


		o_vort = laplace(o_vort, r)*r;
		o_div = laplace(o_div, r)*r;
	}


	/**
	 * Convert vorticity/divergence field to u,v velocity field
	 */
	void robert_vortdiv_to_uv(
			const SphereData_SpectralComplex &i_vrt,
			const SphereData_SpectralComplex &i_div,
			SphereData_PhysicalComplex &o_u,
			SphereData_PhysicalComplex &o_v,

			double i_radius = -1

	)	const
	{
		double ir;

		if (i_radius == -1)
			ir = this->ir;
		else
			ir = 1.0/i_radius;

		i_vrt.request_data_spectral();
		i_div.request_data_spectral();

		SphereData_SpectralComplex psi = inv_laplace(i_vrt, i_radius)*ir;
		SphereData_SpectralComplex chi = inv_laplace(i_div, i_radius)*ir;

		psi.request_data_spectral();
		chi.request_data_spectral();

		SHsphtor_to_spat_cplx(
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
			const SphereData_SpectralComplex &i_phi,
			SphereData_PhysicalComplex &o_u,
			SphereData_PhysicalComplex &o_v,

			double i_radius = -1

	)	const
	{
		double ir;

		if (i_radius == -1)
			ir = this->ir;
		else
			ir = 1.0/i_radius;

		SphereData_SpectralComplex psi(sphereDataConfig);
		psi.spectral_set_zero();

		SphereData_SpectralComplex chi = i_phi;
		chi.request_data_spectral();

		SHsphtor_to_spat_cplx(
						sphereDataConfig->shtns,
						psi.spectral_space_data,
						chi.spectral_space_data,
						o_u.physical_space_data,
						o_v.physical_space_data
				);

		o_u *= ir;
		o_v *= ir;
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

};






#endif /* SPHOPERATORS_HPP_ */
