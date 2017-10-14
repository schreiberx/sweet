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

	/*
	 * Compute Arakawa Jacobian
	 *
	 * J(a,b) = da/dx db/dy - da/dy db/dx
	 */
	PlaneData J(
			const PlaneData &a,
			const PlaneData &b
	)
	{
		return op.diff_c_x(a)*op.diff_c_y(b) - op.diff_c_y(a)*op.diff_c_x(b);
	}


	/*
	 * Compute time derivative of Arakawa Jacobian
	 *
	 * J(a,b)_t = (da/dx db/dy - da/dy db/dx)_t
	 */
	PlaneData J_t(
			PlaneData &a,
			PlaneData &b,
			PlaneData &a_t,
			PlaneData &b_t
	)
	{
		return	  op.diff_c_x(a_t)*op.diff_c_y(b)
				+ op.diff_c_x(a)*op.diff_c_y(b_t)
				- op.diff_c_y(a_t)*op.diff_c_x(b)
				- op.diff_c_y(a)*op.diff_c_x(b_t);
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
		 * Overwrite domain size
		 * Page 179, right column
		 */
		simVars.sim.domain_size[0] = 2.0*M_PI*k0;
		simVars.sim.domain_size[1] = simVars.sim.domain_size[0];

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
		 * Prepare other values
		 */
		// Rossby number
		double R = simVars.swe_polvani.r;
		// Froude number
		double F = simVars.swe_polvani.f;

		double B = (R*R)/(F*F);

#if 0
		std::cout << R << std::endl;
		std::cout << F << std::endl;
		std::cout << B << std::endl;
#endif

		PlaneData laplace = op.diff2_c_x + op.diff2_c_y;

		/*
		 * Compute height
		 */

		PlaneData lap_h = op.diff2(psi) + 2.0*R*J(op.diff_c_x(psi), op.diff_c_y(psi));
		PlaneData h = lap_h.spectral_div_element_wise(laplace);

		o_h = h;
		o_u.spectral_set_zero();
		o_v.spectral_set_zero();

	}
};

#endif
