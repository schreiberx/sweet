/*
 * SWEPolvani.hpp
 *
 *  Created on: 14 Oct 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SWE_PLANE_POLVANI_HPP_
#define SWE_PLANE_POLVANI_HPP_


#include <stdlib.h>
#include <cmath>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>

#if SWEET_THREADING_SPACE
	#include <omp.h>
#endif


/**
 * Implement initial conditions of Polvani benchmark
 * "Coherent structures of shallow-water turbulence"
 */
class SWE_bench_Polvani
{
	SimulationVariables &simVars;

	PlaneOperators &op;



	/*
	 * Start with default parameters
	 */

	// Page 179, left column
	double k0 = 14.0;
	double m = 25;

	double normalize_ek = 1.0;

	double Ek(int ka[2])
	{
		double k = std::sqrt((double)ka[0]*(double)ka[0] + (double)ka[1]*(double)ka[1]);
		return std::pow(k, m/2)/std::pow(k+k0, m);
	}


	void setup_stream(
			PlaneData_Spectral &o_psi
	)
	{
#if SWEET_THREADING_SPACE
		/*
		 * Avoid race conditions in random number generator!
		 */
		int max_threads;

#pragma omp parallel
#pragma omp master
		max_threads = omp_get_num_threads();

		omp_set_num_threads(1);
#endif

		/*
		 * We normalize the energy spectrum to a maximum value of 1
		 */
//		double max_energy_amplitude = 0;

		/*
		 * Determine scaling factor to assure that max of E_k = 1
		 *
		 * http://www.wolframalpha.com/input/?i=solve+diff(k%5E(m%2F2)%2F(k%2Bk_0)%5Em,k)%3D0
		 *
		 * => Maximum at k_0
		 */
		double max_ek = std::pow(k0, m/2.0)/std::pow(2.0*k0, m);

		double scale_ek = 1.0/max_ek;
		std::cout << "POLVANI: Using scaling factor for ek of " << scale_ek << std::endl;

		double scale = o_psi.planeDataConfig->spectral_data_size[0]*o_psi.planeDataConfig->spectral_data_size[1];

		//o_psi.spectral_update_lambda_modes(
		o_psi.spectral_update_lambda(
			[&](int k0, int k1, std::complex<double> &o_data)
			{
#if SWEET_DEBUG && SWEET_THREADING_SPACE
#pragma omp parallel
#pragma omp master
				{
					if (omp_get_num_threads() > 1)
						SWEETError("THREADING MUST BE DEACTIVATED HERE BECAUSE OF RACE CONDITIONS!");
				}

#endif
				int ka[2] = {k0, k1};


//				double energy_amplitude = Ek(ka)*scale_ek;7
				double energy_amplitude = Ek(ka)*normalize_ek;

				double phase = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

//				max_energy_amplitude = std::max(energy_amplitude, max_energy_amplitude);

				/*
				 * Rescale amplitude
				 *
				 * 1/2 * k^2 |\psi_k|^2 = Ek(k)
				 *
				 * \psi_k = SQRT( Ek(k) * 2 / (k^2))
				 */
				double k = std::sqrt((double)k0*(double)k0 + (double)k1*(double)k1);

				if (k == 0.0)
					energy_amplitude = 0;
				else
					energy_amplitude = std::sqrt(energy_amplitude*2.0 / (k*k));

				/*
				 * https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
				 *
				 * Initialize with sin/cos.
				 */
				o_data = {std::sin(phase), std::cos(phase)};
				o_data *= energy_amplitude;

				/*
				 * spectral space scaling factor
				 */
				o_data *= scale;
			}
		);

//		o_psi *= 1.0/max_energy_amplitude;

#if SWEET_THREADING_SPACE
	omp_set_num_threads(max_threads);
#endif
	}


public:
	SWE_bench_Polvani(
		SimulationVariables &io_simVars,
		PlaneOperators &io_op
	)	:
		simVars(io_simVars),
		op(io_op)
	{
	}

	double R;
	double B;
	double F;

	void setup(
			PlaneData_Spectral &o_h,
			PlaneData_Spectral &o_u,
			PlaneData_Spectral &o_v
	)
	{
		/*
		 * Prepare other values
		 */
		// Rossby number
		R = simVars.swe_polvani.r;

		// Froude number
		F = simVars.swe_polvani.f;

		// Burger number
		// Equation (2.2)
		B = (R*R)/(F*F);

#if 0
		std::cout << R << std::endl;
		std::cout << F << std::endl;
		std::cout << B << std::endl;
#endif


		/*
		 * Overwrite domain size
		 * Page 179, right column
		 */
		simVars.sim.plane_domain_size[0] = 2.0*M_PI*k0;
		simVars.sim.plane_domain_size[1] = simVars.sim.plane_domain_size[0];


		/*
		 * Equation (2.3.a)
		 * => Infer f0 and gravitation
		 */
		simVars.sim.plane_rotating_f0 = 1.0/R;
		//simVars.sim.gravitation = 1.0/R;
		simVars.sim.gravitation = 1.0/(F*F);

		//simVars.sim.h0 = 1.0/R*B;
		simVars.sim.h0 = 1.0;

		simVars.outputConfig();

		std::cout << "******************* WARNING ***********************" << std::endl;
		std::cout << "POLVANI Benchmark setup:" << std::endl;
		std::cout << "Updated domain size!" << std::endl;
		std::cout << "******************* WARNING ***********************" << std::endl;

		/*
		 * update domain size
		 */
		op.setup(simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);


//		normalize_ek = 1.5967e+21;
		normalize_ek = 1.59373326082e+21;

		double target_rms = 1.0/sqrt(2.0);

		/*
		 * Use a shitty zero-crossing finder here
		 *
		 * TODO: Use Newton solver
		 */
		for (int k = 0; k < 1000; k++)
		{
			/*
			 * Manually reset the seed
			 */
			if (simVars.benchmark.random_seed >= 0)
				srandom(simVars.benchmark.random_seed);

			std::cout << "POLVANI"<< std::endl;
			std::cout << "POLVANI ITERATION " << k << std::endl;
			std::cout << "POLVANI NORMALIZE_EK = " << normalize_ek << std::endl;

			/*
			 * We first compute the slope numerically
			 */

			/*
			 * This is the distance between two sampling points to determine the slope
			 */
			double d_rms = normalize_ek * 1e-8;

			/*
			 * First sample at x_0
			 */
			double x_0 = normalize_ek;
			double rms_0;
			{
				setup_inner_iter(o_h, o_u, o_v);
				double rms_u = o_u.toPhys().physical_reduce_rms();
				double rms_v = o_v.toPhys().physical_reduce_rms();
				rms_0 = 0.5*(rms_u + rms_v);
			}

			std::cout << "POLVANI x0 = " << x_0 << "\tRMS = " << rms_0 << std::endl;

			/*
			 * Second sample at x_1
			 */
			double x_1 = x_0+d_rms;
			normalize_ek = x_1;
			double rms_1;
			{
				setup_inner_iter(o_h, o_u, o_v);
				double rms_u = o_u.toPhys().physical_reduce_rms();
				double rms_v = o_v.toPhys().physical_reduce_rms();
				rms_1 = 0.5*(rms_u + rms_v);
			}

			std::cout << "POLVANI x1 = " << x_1 << "\tRMS = " << rms_1 << std::endl;

			double diff = (rms_1 - rms_0) / d_rms;

			std::cout << "POLVANI dy/dx = " << diff << std::endl;

			double omega = 100.0;

			normalize_ek = normalize_ek - omega*(rms_0 - target_rms)/diff;

			double error = std::abs(rms_0 - target_rms);
			std::cout << "POLVANI " << k << ": " << rms_0 << "\tERROR: " << error << std::endl;
			std::cout << std::endl;

			if (error < 1e-2)
				break;
		}
	}



	void setup_inner_iter(
			PlaneData_Spectral &o_h,
			PlaneData_Spectral &o_u,
			PlaneData_Spectral &o_v
	)
	{

		PlaneData_Spectral psi(o_h.planeDataConfig);

		/*
		 * Prepare laplace operator
		 */
		PlaneData_Spectral laplace = op.diff2_c_x + op.diff2_c_y;

		/*
		 * Setup stream function
		 * Equation (2.6)
		 */
		setup_stream(psi);

//		psi.print_physicalArrayData();
//		std::cout << lap_h.reduce_maxAbs() << std::endl;

		/*
		 * Compute height
		 * Solve equation (2.5b)
		 */

		PlaneData_Spectral lap_h = op.diff2(psi) + 2.0*R*op.J(op.diff_c_x(psi), op.diff_c_y(psi));
		PlaneData_Spectral h = lap_h.spectral_div_element_wise(laplace);

		/*
		 * Setup chi
		 */
		PlaneData_Spectral chi(o_h.planeDataConfig);
		chi.spectral_set_zero();

		PlaneData_Spectral psi_t(o_h.planeDataConfig);

		/*
		 * Iteratively solve for chi
		 */
		double R_1 = 1.0/R;
		double diff;
		for (int i = 0; i < 50; i++)
		{
			/*
			 * Solve equation (2.5a)
			 */
			PlaneData_Spectral laplace_psi_t =
					  op.J(psi, laplace(psi))
					+ R_1*(laplace(chi))
					+ op.div(	laplace(chi)*op.diff_c_x(chi),
								laplace(chi)*op.diff_c_y(chi)
						)
					;

			PlaneData_Spectral psi_t = laplace_psi_t.spectral_div_element_wise(laplace);

			/*
			 * Solve equation (2.5c)
			 */
			PlaneData_Spectral stuff_chi =
					- op.J(psi, laplace(chi))
					+ laplace(op.J(psi, h))
					+ 2.0*R*op.J_t(psi, psi, op.diff_c_x(psi), op.diff_c_y(psi))
					- op.div(	laplace(chi)*op.diff_c_x(chi),
								laplace(chi)*op.diff_c_y(chi)
						)
					+ laplace(op.diff_c_x(h*op.diff_c_x(chi)) + op.diff_c_y(h*op.diff_c_y(chi)));

			/*
			 * Setup identity operator for Helmholtz solver
			 */
			PlaneData_Spectral I(o_h.planeDataConfig);
			I.spectral_set_zero();
			I.spectral_addScalarAll(1.0);

			PlaneData_Spectral lhs = R_1*(I - laplace*B);
			PlaneData_Spectral new_chi = stuff_chi.spectral_div_element_wise(laplace).spectral_div_element_wise(lhs);

			diff = (new_chi-chi).spectral_reduce_max_abs();
			std::cout << i << ": chi update = " << diff << std::endl;

			chi = new_chi;

			if (diff < 1e-10)
				break;
		}

		if (diff >= 1e-10)
		{
			SWEETError("No convergence for Polvani initial conditions reached");
		}


		/*
		 * Convert height field
		 * See page 178, h* = H(1+RB^-1 h)
		 */
		// total height
		h = simVars.sim.h0*(1.0+R/B*h);

		// perturbation
		h = h-simVars.sim.h0;

		o_h = h;


		/*
		 * Compute velocity
		 */
		double eps = R*std::max(1.0, R)/std::max(1.0, B);

		// Equation (2.4)
		o_u = -op.diff_c_y(psi) + eps*op.diff_c_x(chi);
		o_v = op.diff_c_x(psi) + eps*op.diff_c_y(chi);

		double chi_rms = chi.toPhys().physical_reduce_rms();
		double psi_rms = psi.toPhys().physical_reduce_rms();

		std::cout << "POLVANI: chi_rms = " << chi_rms << std::endl;
		std::cout << "POLVANI: psi_rms = " << psi_rms << std::endl;
		double chi_rms_psi_rms = chi_rms / psi_rms;
		std::cout << "POLVANI: chi_rms / psi_rms = " << chi_rms_psi_rms << std::endl;
	}
};








#endif
