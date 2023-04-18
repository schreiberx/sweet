/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "SWE_Plane_Mori_Zwanzig_TS_l_direct.hpp"

#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/plane/PlaneOperatorsComplex.hpp>

#include <sweet/core/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>
#include <sweet/core/plane/PlaneStaggering.hpp>



bool SWE_Plane_Mori_Zwanzig_TS_l_direct::setup(
		sweet::PlaneOperators *io_ops,
		const std::string &i_function_name
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	assert(shackPlaneDataOps != nullptr);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		planeDataGridMapping.setup(shackPlaneDataOps, ops->planeDataConfig);

	expFunctions.setup(i_function_name);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(expFunctions);

	this->normal_modes.setup();

	return true;
}


bool SWE_Plane_Mori_Zwanzig_TS_l_direct::setup(
		sweet::PlaneOperators *io_ops
)
{
	return setup(io_ops, "phi0");
}


void SWE_Plane_Mori_Zwanzig_TS_l_direct::runTimestep(
		sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (shackPlaneDataOps->space_grid_use_c_staggering)
		run_timestep_cgrid(io_h_pert, io_u, io_v, i_dt, i_simulation_timestamp);
	else
		run_timestep_agrid(io_h_pert, io_u, io_v, i_dt, i_simulation_timestamp);
}


/**
 * Computation of analytical solution on staggered grid
 */
void SWE_Plane_Mori_Zwanzig_TS_l_direct::run_timestep_cgrid(
		sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,		///< prognostic variables
		sweet::PlaneData_Spectral &io_v,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	// For output, variables need to be on unstaggered A-grid
	sweet::PlaneData_Physical t_u(io_h_pert.planeDataConfig);
	sweet::PlaneData_Physical t_v(io_h_pert.planeDataConfig);

	if (!shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Expected staggering");

	planeDataGridMapping.mapCtoA_u(io_u.toPhys(), t_u);
	planeDataGridMapping.mapCtoA_v(io_v.toPhys(), t_v);

	shackPlaneDataOps->space_grid_use_c_staggering = false;

	sweet::PlaneData_Spectral t_u_spec(io_h_pert.planeDataConfig);
	sweet::PlaneData_Spectral t_v_spec(io_h_pert.planeDataConfig);
	t_u_spec.loadPlaneDataPhysical(t_u);
	t_v_spec.loadPlaneDataPhysical(t_v);

	run_timestep_agrid(
			io_h_pert, t_u_spec, t_v_spec,
			i_dt, i_simulation_timestamp
	);

	shackPlaneDataOps->space_grid_use_c_staggering = true;

	planeDataGridMapping.mapAtoC_u(t_u_spec.toPhys(), t_u);
	planeDataGridMapping.mapAtoC_v(t_v_spec.toPhys(), t_v);

	io_u.loadPlaneDataPhysical(t_u);
	io_v.loadPlaneDataPhysical(t_v);
}



void SWE_Plane_Mori_Zwanzig_TS_l_direct::run_timestep_agrid(
		sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{

	run_timestep_agrid_planedata(io_h_pert, io_u, io_v, i_dt, i_simulation_timestamp);
}


#if SWEET_USE_PLANE_SPECTRAL_SPACE

/**
 * This method computes the analytical solution based on the given initial values.
 *
 * See Embid/Madja/1996, Terry/Beth/2014, page 16
 * and
 * 		doc/swe_solution_for_L/sympy_L_spec_decomposition.py
 * for the dimension full formulation.
 */
void SWE_Plane_Mori_Zwanzig_TS_l_direct::run_timestep_agrid_planedata(
		sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported");


	typedef std::complex<T> complex;
	complex I(0.0, 1.0);

	T dt = i_dt;


	// compute Q*phin(Dt*Lambda)*Q^{-1}
#if SWEET_THREADING_SPACE
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
#endif
	for (std::size_t ik1 = 0; ik1 < io_h_pert.planeDataConfig->spectral_data_size[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < io_h_pert.planeDataConfig->spectral_data_size[0]; ik0++)
		{

			complex eigenvalues[3];
			complex eigenvectors[3][3];

			T k1;
			if (ik1 < io_h_pert.planeDataConfig->spectral_data_size[1]/2)
				k1 = (T)ik1;
			else
				k1 = (T)((int)ik1-(int)io_h_pert.planeDataConfig->spectral_data_size[1]);

			T k0 = (T)ik0;


			normal_modes.eigendecomposition(k0, k1, eigenvalues, eigenvectors);

			/*
			 * Invert Eigenvalue matrix
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

			// Compute Q * phin(Dt * Lambda)
			complex v_lambda[3][3];
			for (int i = 0; i < 3; i++)
			{

				std::complex<T> &lam = eigenvalues[i];

				std::complex<T> K = expFunctions.eval(lam*dt);
				for (int j = 0; j < 3; j++)
					v_lambda[j][i] = eigenvectors[j][i] * K;
			}

			complex U[3];
			U[0] = io_h_pert.spectral_get(ik1, ik0);
			U[1] = io_u.spectral_get(ik1, ik0);
			U[2] = io_v.spectral_get(ik1, ik0);

			complex U_copy[3];
			for (int k = 0; k < 3; k++)
			{
				U_copy[k] = U[k];
				U[k] = 0.0;
			}


			// Compute  [Q*phin(Dt*Lambda)]*Q^{-1}
			for (int j = 0; j < 3; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					std::complex<double> d = 0.;
					for (int k = 0; k < 3; k++)
						d += v_lambda[j][k] * eigenvectors_inv[k][i];

					U[j] += d * U_copy[j];
				}
			}

#if SWEET_QUADMATH
			std::complex<double> tmp0(U[0].real(), U[0].imag());
			io_h_pert.spectral_set(ik1, ik0, tmp0);

			std::complex<double> tmp1(U[1].real(), U[1].imag());
			io_u.spectral_set(ik1, ik0, tmp1);

			std::complex<double> tmp2(U[2].real(), U[2].imag());
			io_v.spectral_set(ik1, ik0, tmp2);
#else
			io_h_pert.spectral_set(ik1, ik0, U[0]);
			io_u.spectral_set(ik1, ik0, U[1]);
			io_v.spectral_set(ik1, ik0, U[2]);
#endif
		}
	}
	
	io_h_pert.spectral_zeroAliasingModes();
	io_u.spectral_zeroAliasingModes();
	io_v.spectral_zeroAliasingModes();
}

#endif


