/*
 * SWE_Plane_TS_l_direct.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Plane_TS_l_direct.hpp"
#include <sweet/plane/PlaneDataComplex.hpp>

#include <sweet/plane/PlaneOperatorsComplex.hpp>

#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>

/**
 * This method computes the analytical solution based on the given initial values.
 *
 * See Embid/Madja/1996, Terry/Beth/2014, page 16
 * and
 * 		doc/swe_solution_for_L/sympy_L_spec_decomposition.py
 * for the dimensionful formulation.
 *
 * Don't use this function too frequently, since it always computes
 * the required coefficients on-the-fly which is expensive.
 */
void SWE_Plane_TS_l_direct::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt < 0)
		FatalError("SWE_Plane_TS_l_direct: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time - i_simulation_timestamp;



	typedef std::complex<double> complex;

	double eta_bar = simVars.sim.h0;
	double g = simVars.sim.gravitation;
	double f = simVars.sim.f0;
	complex I(0.0,1.0);

	PlaneDataComplex i_h = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_h);
	PlaneDataComplex i_u = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_u);
	PlaneDataComplex i_v = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_v);

	PlaneDataComplex o_h(io_h.planeDataConfig);
	PlaneDataComplex o_u(io_h.planeDataConfig);
	PlaneDataComplex o_v(io_h.planeDataConfig);

	double s0 = simVars.sim.domain_size[0];
	double s1 = simVars.sim.domain_size[1];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	o_h.spectral_space_data_valid = true;
	o_h.physical_space_data_valid = false;

	o_u.spectral_space_data_valid = true;
	o_u.physical_space_data_valid = false;

	o_v.spectral_space_data_valid = true;
	o_v.physical_space_data_valid = false;
#endif

	for (std::size_t ik1 = 0; ik1 < i_h.planeDataConfig->spectral_complex_data_size[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < i_h.planeDataConfig->spectral_complex_data_size[0]; ik0++)
		{
			if (ik0 == i_h.planeDataConfig->spectral_complex_data_size[0]/2 || ik1 == i_h.planeDataConfig->spectral_complex_data_size[1]/2)
			{
				o_h.p_spectral_set(ik1, ik0, 0, 0);
				o_u.p_spectral_set(ik1, ik0, 0, 0);
				o_v.p_spectral_set(ik1, ik0, 0, 0);
			}

			complex U_hat[3];
			U_hat[0] = i_h.spectral_get(ik1, ik0);
			U_hat[1] = i_u.spectral_get(ik1, ik0);
			U_hat[2] = i_v.spectral_get(ik1, ik0);

			double k0, k1;
			if (ik0 < i_h.planeDataConfig->spectral_complex_data_size[0]/2)
				k0 = (double)ik0;
			else
				k0 = (double)((int)ik0-(int)i_h.planeDataConfig->spectral_complex_data_size[0]);

			if (ik1 < i_h.planeDataConfig->spectral_complex_data_size[1]/2)
				k1 = (double)ik1;
			else
				k1 = (double)((int)ik1-(int)i_h.planeDataConfig->spectral_complex_data_size[1]);

			/*
			 * dimensionful formulation
			 * see doc/swe_solution_for_L
			 */

			double H0 = eta_bar;

			//////////////////////////////////////
			// GENERATED CODE START
			//////////////////////////////////////
			complex eigenvalues[3];
			complex eigenvectors[3][3];

			if (k0 == 0 && k1 == 0)
			{
				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*f;
				eigenvalues[2] = f;

				eigenvectors[0][0] = 1.00000000000000;
				eigenvectors[0][1] = 0.0;
				eigenvectors[0][2] = 0.0;
				eigenvectors[1][0] = 0.0;
				eigenvectors[1][1] = -1.0*I;
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 0.0;
				eigenvectors[2][1] = I;
				eigenvectors[2][2] = 1.00000000000000;
			}
			else if (k0 == 0)
			{
//					complex wg = std::sqrt((complex)s0*s0*(f*f*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1));

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvalues[2] = -1.0*I*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k1*k1 - 1.0*f*f*s1*s1);

				eigenvectors[0][0] = (1.0/2.0)*I*1.0/M_PI*f*1.0/g*1.0/k1*s1;
				eigenvectors[0][1] = 1.00000000000000;
				eigenvectors[0][2] = 0.0;
				eigenvectors[1][0] = -2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][1] = -1.0*I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][1] = I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][2] = 1.00000000000000;
			}
			else if (k1 == 0)
			{
//					complex wg = std::sqrt((complex)s1*s1*(f*f*s0*s0 + 4.0*M_PI*M_PI*g*g*k0*k0));

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
				eigenvalues[2] = -1.0*I*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);

				eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
				eigenvectors[0][1] = 0.0;
				eigenvectors[0][2] = 1.00000000000000;
				eigenvectors[1][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
				eigenvectors[1][1] = -1.0*I*1.0/f*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
				eigenvectors[2][1] = 1.0/f*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);
				eigenvectors[2][2] = 1.00000000000000;
			}
			else
			{
//					complex K2 = M_PI*M_PI*k0*k0 + M_PI*M_PI*k1*k1;
				complex w = std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);

//					complex wg = std::sqrt((complex)f*f*s0*s0*s1*s1 + 4.0*M_PI*M_PI*g*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1*s0*s0);
				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);
				eigenvalues[2] = -1.0*I*1.0/s0*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 - 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 - 1.0*f*f*s0*s0*s1*s1);

				eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
				eigenvectors[0][1] = -1.0*1.0/k0*k1*s0*1.0/s1;
				eigenvectors[0][2] = 1.00000000000000;
				eigenvectors[1][0] = 2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 + f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 - 1.0*I*f*w);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = -2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 - 1.0*f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 + I*f*w);
				eigenvectors[2][2] = 1.00000000000000;
			}




			//////////////////////////////////////
			// GENERATED CODE END
			//////////////////////////////////////


			if (f == 0)
			{
				/*
				 * override if f == 0, see ./sympy_L_spec_decomposition.py executed with LNr=4
				 */
				if (k0 != 0 || k1 != 0)
				{
					double K2 = K2;

					eigenvalues[0] = 0.0;
					eigenvalues[1] = -2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);
					eigenvalues[2] = 2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);

					eigenvectors[0][0] = 0.0;
					eigenvectors[0][1] = -1.0*k1/sqrt(k0*k0 + k1*k1);
					eigenvectors[0][2] = k0/sqrt(k0*k0 + k1*k1);
					eigenvectors[1][0] = -1.0*sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[1][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[1][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][0] = sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
				}
				else
				{

					eigenvalues[0] = 0.0;
					eigenvalues[1] = 0.0;
					eigenvalues[2] = 0.0;

					eigenvectors[0][0] = 1.00000000000000;
					eigenvectors[0][1] = 0.0;
					eigenvectors[0][2] = 0.0;
					eigenvectors[1][0] = 0.0;
					eigenvectors[1][1] = 1.00000000000000;
					eigenvectors[1][2] = 0.0;
					eigenvectors[2][0] = 0.0;
					eigenvectors[2][1] = 0.0;
					eigenvectors[2][2] = 1.00000000000000;
				}
			}


			/*
			 * Compute inverse of Eigenvectors.
			 * This generalizes to the case that the Eigenvectors are not orthonormal.
			 */
			complex eigenvectors_inv[3][3];

			eigenvectors_inv[0][0] =  (eigenvectors[1][1]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][1]);
			eigenvectors_inv[0][1] = -(eigenvectors[0][1]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][1]);
			eigenvectors_inv[0][2] =  (eigenvectors[0][1]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][1]);

			eigenvectors_inv[1][0] = -(eigenvectors[1][0]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][0]);
			eigenvectors_inv[1][1] =  (eigenvectors[0][0]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][0]);
			eigenvectors_inv[1][2] = -(eigenvectors[0][0]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][0]);

			eigenvectors_inv[2][0] =  (eigenvectors[1][0]*eigenvectors[2][1] - eigenvectors[1][1]*eigenvectors[2][0]);
			eigenvectors_inv[2][1] = -(eigenvectors[0][0]*eigenvectors[2][1] - eigenvectors[0][1]*eigenvectors[2][0]);
			eigenvectors_inv[2][2] =  (eigenvectors[0][0]*eigenvectors[1][1] - eigenvectors[0][1]*eigenvectors[1][0]);

			complex s = eigenvectors[0][0]*eigenvectors_inv[0][0] + eigenvectors[0][1]*eigenvectors_inv[1][0] + eigenvectors[0][2]*eigenvectors_inv[2][0];

			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++)
					eigenvectors_inv[j][i] /= s;


			// check
			for (int j = 0; j < 3; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					if (
							std::isnan(eigenvectors[j][i].real()) || std::isinf(eigenvectors[j][i].real()) != 0	||
							std::isnan(eigenvectors[j][i].imag()) || std::isinf(eigenvectors[j][i].imag()) != 0
					)
					{
						std::cerr << "Invalid number in Eigenvector " << j << " detected: " << eigenvectors[j][0] << ", " << eigenvectors[j][1] << ", " << eigenvectors[j][2] << std::endl;
					}

					if (
							std::isnan(eigenvectors_inv[j][i].real()) || std::isinf(eigenvectors_inv[j][i].real()) != 0	||
							std::isnan(eigenvectors_inv[j][i].imag()) || std::isinf(eigenvectors_inv[j][i].imag()) != 0
					)
					{
						std::cerr << "Invalid number in inverse of Eigenvector " << j << " detected: " << eigenvectors_inv[j][0] << ", " << eigenvectors_inv[j][1] << ", " << eigenvectors_inv[j][2] << std::endl;
					}
				}
			}

			/*
			 * Solve based on previously computed data.
			 * Note, that this data can be also precomputed and reused every time.
			 */
			complex UEV0_sp[3];
			for (int k = 0; k < 3; k++)
			{
				UEV0_sp[k] = {0, 0};
				for (int j = 0; j < 3; j++)
					UEV0_sp[k] += eigenvectors_inv[j][k] * U_hat[j];
			}

			complex omega[3];
			omega[0] = std::exp(-I*eigenvalues[0]*i_fixed_dt);
			omega[1] = std::exp(-I*eigenvalues[1]*i_fixed_dt);
			omega[2] = std::exp(-I*eigenvalues[2]*i_fixed_dt);

			complex U_hat_sp[3];
			for (int k = 0; k < 3; k++)
			{
				U_hat_sp[k] = {0, 0};
				for (int j = 0; j < 3; j++)
					U_hat_sp[k] += eigenvectors[j][k] * omega[j] * UEV0_sp[j];
			}

			o_h.p_spectral_set(ik1, ik0, U_hat_sp[0]);
			o_u.p_spectral_set(ik1, ik0, U_hat_sp[1]);
			o_v.p_spectral_set(ik1, ik0, U_hat_sp[2]);
		}
	}

	io_h = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_h);
	io_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_u);
	io_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_v);

	o_dt = i_fixed_dt;
}





SWE_Plane_TS_l_direct::SWE_Plane_TS_l_direct(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
}



SWE_Plane_TS_l_direct::~SWE_Plane_TS_l_direct()
{
}

