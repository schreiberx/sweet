/*
 * SWEPolvani.hpp
 *
 *  Created on: 14 Oct 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SWE_PLANE_POLVANI_HPP_
#define SWE_PLANE_POLVANI_HPP_


#include <stdlib.h>
#include <sweet/sweetmath.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>

#if SWEET_SPACE_THREADING
	#include <omp.h>
#endif


/**
 * Implement initial conditions of Polvani benchmark
 * "Coherent structures of shallow-water turbulence"
 */
class SWEPolvani
{
	SimulationVariables &simVars;

	PlaneOperators &op;



	/*
	 * Start with default parameters
	 */

	// Page 179, left column
	double k0 = 14.0;
	double m = 25;


	double Ek(int ka[2])
	{
		double k = std::abs(ka[0]) + std::abs(ka[1]);
		return std::pow(k, m/2)/std::pow(k+k0, m);
	}


	void setup_stream(
			PlaneData &o_psi
	)
	{
#if SWEET_SPACE_THREADING
		/*
		 * Avoid race conditions in random number generator!
		 */
		int max_threads = omp_get_max_threads();
		omp_set_num_threads(1);
#endif

		/*
		 * We normalize the energy spectrum to a maximum value of 1
		 */
		double max_amplitude = 0;

		o_psi.spectral_update_lambda_modes(
			[&](int k0, int k1, std::complex<double> &o_data)
			{
#if SWEET_DEBUG && SWEET_SPACE_THREADING
				if (omp_get_num_threads() > 1)
					FatalError("THREADING MUST BE DEACTIVATED HERE BECAUSE OF RACE CONDITIONS!");
#endif
				int ka[2] = {k0, k1};

				double k = std::abs(k0) + std::abs(k1);

				double amplitude = Ek(ka);
				double phase = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

				/*
				 * Rescale amplitude
				 *
				 * 1/2 * k^2 |\psi_k|^2 = Ek(k)
				 *
				 * \psi_k = SQRT( Ek(k) * 2 / (k^2))
				 */
				if (k == 0)
					amplitude = 0;
				else
					amplitude = std::sqrt(amplitude*2.0 / (k*k));

				max_amplitude = std::max(amplitude, max_amplitude);

				/*
				 * https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
				 *
				 * Initialize with sin/cos.
				 */
				o_data = {std::sin(phase), std::cos(phase)};
				o_data *= amplitude;

				// Not necessary, because of scaling to max amplitude
//				o_data *= o_psi.planeDataConfig->physical_array_data_number_of_elements;
			}
		);

		o_psi *= 1.0/max_amplitude;

#if SWEET_SPACE_THREADING
	omp_set_num_threads(max_threads);
#endif
	}


public:
	SWEPolvani(
		SimulationVariables &io_simVars,
		PlaneOperators &io_op
	)	:
		simVars(io_simVars),
		op(io_op)
	{
	}



	void setup(
			PlaneData &o_h,
			PlaneData &o_u,
			PlaneData &o_v
	)
	{
		/*
		 * Prepare other values
		 */
		// Rossby number
		double R = simVars.swe_polvani.r;
		// Froude number
		double F = simVars.swe_polvani.f;

		// Burger number
		// Equation (2.2)
		double B = (R*R)/(F*F);

#if 0
		std::cout << R << std::endl;
		std::cout << F << std::endl;
		std::cout << B << std::endl;
#endif


		/*
		 * Overwrite domain size
		 * Page 179, right column
		 */
		simVars.sim.domain_size[0] = 2.0*M_PI*k0;
		simVars.sim.domain_size[1] = simVars.sim.domain_size[0];


		/*
		 * Equation (2.3.a)
		 * => Infer f0 and gravitation
		 */
		simVars.sim.f0 = 1.0/R;
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
		op.setup(simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

		PlaneData psi(o_h.planeDataConfig);

		/*
		 * Setup stream function
		 * Equation (2.6)
		 */
		setup_stream(psi);

//		psi.print_physicalArrayData();
//		std::cout << lap_h.reduce_maxAbs() << std::endl;

		/*
		 * Prepare laplace operator
		 */
		PlaneData laplace = op.diff2_c_x + op.diff2_c_y;

		/*
		 * Compute height
		 * Solve equation (2.5b)
		 */

		PlaneData lap_h = op.diff2(psi) + 2.0*R*op.J(op.diff_c_x(psi), op.diff_c_y(psi));
		PlaneData h = lap_h.spectral_div_element_wise(laplace);

		/*
		 * Setup chi
		 */
		PlaneData chi(o_h.planeDataConfig);
		chi.spectral_set_zero();

		PlaneData psi_t(o_h.planeDataConfig);

		/*
		 * Iteratively solve for chi
		 */
		double R_1 = 1.0/R;
		double diff;
		for (int i = 0; i < 10; i++)
		{
			/*
			 * Solve equation (2.5a)
			 */
			PlaneData laplace_psi_t =
					  op.J(psi, laplace(psi))
					+ R_1*(laplace(chi))
					+ op.div(	laplace(chi)*op.diff_c_x(chi),
								laplace(chi)*op.diff_c_y(chi)
						)
					;

			PlaneData psi_t = laplace_psi_t.spectral_div_element_wise(laplace);

			/*
			 * Solve equation (2.5c)
			 */
			PlaneData stuff_chi =
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
			PlaneData I(o_h.planeDataConfig);
			I.spectral_set_zero();
			I.spectral_addScalarAll(1.0);

			PlaneData lhs = R_1*(I - laplace*B);
			PlaneData new_chi = stuff_chi.spectral_div_element_wise(laplace).spectral_div_element_wise(lhs);

			diff = (new_chi-chi).reduce_maxAbs();
			std::cout << i << ": chi update = " << diff << std::endl;

			chi = new_chi;
		}

		if (diff > 1e-10)
		{
			FatalError("No convergence for Polvani initial conditions reached");
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

	}
};








#endif
