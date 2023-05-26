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


bool SWE_Plane_Mori_Zwanzig_TS_l_direct::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{

	PDESWEPlaneMoriZwanzigTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}


bool SWE_Plane_Mori_Zwanzig_TS_l_direct::setup(
		sweet::PlaneOperators *io_ops,
		const std::string &i_function_name
)
{
	ops = io_ops;

	assert(shackPlaneDataOps != nullptr);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		planeDataGridMapping.setup(shackPlaneDataOps, ops->planeDataConfig);

	expFunctions.setup(i_function_name);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(expFunctions);

	this->normal_modes.setup(
					shackPDESWEPlane->plane_rotating_f0,
					shackPDESWEPlane->h0,
					shackPDESWEPlane->gravitation,
					shackPlaneDataOps->plane_domain_size[0],
					shackPlaneDataOps->plane_domain_size[1],
					ops->planeDataConfig
				);


	return true;
}


bool SWE_Plane_Mori_Zwanzig_TS_l_direct::setup(
		sweet::PlaneOperators *io_ops
)
{
	return setup(io_ops, "phi0");
}


void SWE_Plane_Mori_Zwanzig_TS_l_direct::runTimestep(

		sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables


		double i_dt,
		double i_simulation_timestamp
)
{
	if (shackPlaneDataOps->space_grid_use_c_staggering)
		run_timestep_cgrid(
					io_h_pert_SP, io_u_SP, io_v_SP,
					io_h_pert_SQ, io_u_SQ, io_v_SQ,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt, i_simulation_timestamp
				);
	else
		run_timestep_agrid(
					io_h_pert_SP, io_u_SP, io_v_SP,
					io_h_pert_SQ, io_u_SQ, io_v_SQ,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt, i_simulation_timestamp
				);
}


/**
 * Computation of analytical solution on staggered grid
 */
void SWE_Plane_Mori_Zwanzig_TS_l_direct::run_timestep_cgrid(

		sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	// For output, variables need to be on unstaggered A-grid
	sweet::PlaneData_Physical t_u_SP(io_h_pert_SP.planeDataConfig);
	sweet::PlaneData_Physical t_v_SP(io_h_pert_SP.planeDataConfig);
	sweet::PlaneData_Physical t_u_SQ(io_h_pert_SQ.planeDataConfig);
	sweet::PlaneData_Physical t_v_SQ(io_h_pert_SQ.planeDataConfig);
	sweet::PlaneData_Physical t_u_FQ(io_h_pert_FQ.planeDataConfig);
	sweet::PlaneData_Physical t_v_FQ(io_h_pert_FQ.planeDataConfig);

	if (!shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Expected staggering");

	planeDataGridMapping.mapCtoA_u(io_u_SP.toPhys(), t_u_SP);
	planeDataGridMapping.mapCtoA_v(io_v_SP.toPhys(), t_v_SP);
	planeDataGridMapping.mapCtoA_u(io_u_SQ.toPhys(), t_u_SQ);
	planeDataGridMapping.mapCtoA_v(io_v_SQ.toPhys(), t_v_SQ);
	planeDataGridMapping.mapCtoA_u(io_u_FQ.toPhys(), t_u_FQ);
	planeDataGridMapping.mapCtoA_v(io_v_FQ.toPhys(), t_v_FQ);

	shackPlaneDataOps->space_grid_use_c_staggering = false;

	sweet::PlaneData_Spectral t_u_SP_spec(io_h_pert_SP.planeDataConfig);
	sweet::PlaneData_Spectral t_v_SP_spec(io_h_pert_SP.planeDataConfig);
	sweet::PlaneData_Spectral t_u_SQ_spec(io_h_pert_SQ.planeDataConfig);
	sweet::PlaneData_Spectral t_v_SQ_spec(io_h_pert_SQ.planeDataConfig);
	sweet::PlaneData_Spectral t_u_FQ_spec(io_h_pert_FQ.planeDataConfig);
	sweet::PlaneData_Spectral t_v_FQ_spec(io_h_pert_FQ.planeDataConfig);
	t_u_SP_spec.loadPlaneDataPhysical(t_u_SP);
	t_v_SP_spec.loadPlaneDataPhysical(t_v_SP);
	t_u_SP_spec.loadPlaneDataPhysical(t_u_SQ);
	t_v_SP_spec.loadPlaneDataPhysical(t_v_SQ);
	t_u_SP_spec.loadPlaneDataPhysical(t_u_FQ);
	t_v_SP_spec.loadPlaneDataPhysical(t_v_FQ);

	run_timestep_agrid(
					io_h_pert_SP, t_u_SP_spec, t_v_SP_spec,
					io_h_pert_SQ, t_u_SQ_spec, t_v_SQ_spec,
					io_h_pert_FQ, t_u_FQ_spec, t_v_FQ_spec,
					i_dt, i_simulation_timestamp
			);

	shackPlaneDataOps->space_grid_use_c_staggering = true;

	planeDataGridMapping.mapAtoC_u(t_u_SP_spec.toPhys(), t_u_SP);
	planeDataGridMapping.mapAtoC_v(t_v_SP_spec.toPhys(), t_v_SP);
	planeDataGridMapping.mapAtoC_u(t_u_SP_spec.toPhys(), t_u_SQ);
	planeDataGridMapping.mapAtoC_v(t_v_SP_spec.toPhys(), t_v_SQ);
	planeDataGridMapping.mapAtoC_u(t_u_SP_spec.toPhys(), t_u_FQ);
	planeDataGridMapping.mapAtoC_v(t_v_SP_spec.toPhys(), t_v_FQ);

	io_u_SP.loadPlaneDataPhysical(t_u_SP);
	io_v_SP.loadPlaneDataPhysical(t_v_SP);
	io_u_SQ.loadPlaneDataPhysical(t_u_SQ);
	io_v_SQ.loadPlaneDataPhysical(t_v_SQ);
	io_u_FQ.loadPlaneDataPhysical(t_u_FQ);
	io_v_FQ.loadPlaneDataPhysical(t_v_FQ);
}



void SWE_Plane_Mori_Zwanzig_TS_l_direct::run_timestep_agrid(

		sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{

	run_timestep_agrid_planedata(
					io_h_pert_SP, io_u_SP, io_v_SP,
					io_h_pert_SQ, io_u_SQ, io_v_SQ,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt, i_simulation_timestamp
			);
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

		sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables


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
	for (std::size_t ik1 = 0; ik1 < io_h_pert_SP.planeDataConfig->spectral_data_size[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < io_h_pert_SP.planeDataConfig->spectral_data_size[0]; ik0++)
		{

			complex eigenvalues[3];
			complex eigenvectors[3][3];
			complex eigenvectors_inv[3][3];

			////T k1;
			////if (ik1 < io_h_pert_SP.planeDataConfig->spectral_data_size[1]/2)
			////	k1 = (T)ik1;
			////else
			////	k1 = (T)((int)ik1-(int)io_h_pert_SP.planeDataConfig->spectral_data_size[1]);

			////T k0 = (T)ik0;

			// get eigenvectors matrix, its inverse and eigenvalues
			normal_modes.eigendecomposition(ik0, ik1, eigenvalues, eigenvectors, eigenvectors_inv);

			complex U_SP[3];
			complex U_FQ[3];
			U_SP[0] = io_h_pert_SP.spectral_get(ik1, ik0);
			U_SP[1] = io_u_SP.spectral_get(ik1, ik0);
			U_SP[2] = io_v_SP.spectral_get(ik1, ik0);
			U_FQ[0] = io_h_pert_FQ.spectral_get(ik1, ik0);
			U_FQ[1] = io_u_FQ.spectral_get(ik1, ik0);
			U_FQ[2] = io_v_FQ.spectral_get(ik1, ik0);

			// Compute Qinv * U
			complex UEV_SP[3] = {0., 0., 0.};
			complex UEV_FQ[3] = {0., 0., 0.};
			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
				{
					UEV_SP[k] += eigenvectors_inv[k][j] * U_SP[j];
					UEV_FQ[k] += eigenvectors_inv[k][j] * U_FQ[j];
				}

			// Compute phin(Dt * Lambda) * [Qinv * U]
			for (int k = 0; k < 3; k++)
			{
				std::complex<T> &lam = eigenvalues[k];

				std::complex<T> K = expFunctions.eval(lam*dt);

				UEV_SP[k] = K * UEV_SP[k];
				UEV_FQ[k] = K * UEV_FQ[k];
			}

			// Compute Q * [phin * Qinv * U]
			for (int k = 0; k < 3; k++)
			{
				U_SP[k] = 0.;
				U_FQ[k] = 0.;
			}
			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
				{
					U_SP[k] += eigenvectors[k][j] * UEV_SP[j];
					U_FQ[k] += eigenvectors[k][j] * UEV_FQ[j];
				}

#if SWEET_QUADMATH
			std::complex<double> tmp0(U_SP[0].real(), U_SP[0].imag());
			io_h_pert_SP.spectral_set(ik1, ik0, tmp0);

			std::complex<double> tmp1(U_SP[1].real(), U_SP[1].imag());
			io_u_SP.spectral_set(ik1, ik0, tmp1);

			std::complex<double> tmp2(U_SP[2].real(), U_SP[2].imag());
			io_v_SP.spectral_set(ik1, ik0, tmp2);

			std::complex<double> tmp0(U_FQ[0].real(), U_FQ[0].imag());
			io_h_pert_FQ.spectral_set(ik1, ik0, tmp3);

			std::complex<double> tmp1(U_FQ[1].real(), U_FQ[1].imag());
			io_u_FQ.spectral_set(ik1, ik0, tmp4);

			std::complex<double> tmp2(U_FQ[2].real(), U_FQ[2].imag());
			io_v_FQ.spectral_set(ik1, ik0, tmp5);


#else
			io_h_pert_SP.spectral_set(ik1, ik0, U_SP[0]);
			io_u_SP.spectral_set(ik1, ik0, U_SP[1]);
			io_v_SP.spectral_set(ik1, ik0, U_SP[2]);
			io_h_pert_FQ.spectral_set(ik1, ik0, U_FQ[0]);
			io_u_FQ.spectral_set(ik1, ik0, U_FQ[1]);
			io_v_FQ.spectral_set(ik1, ik0, U_FQ[2]);
#endif
		}
	}

	io_h_pert_SP.spectral_zeroAliasingModes();
	io_u_SP.spectral_zeroAliasingModes();
	io_v_SP.spectral_zeroAliasingModes();
	io_h_pert_FQ.spectral_zeroAliasingModes();
	io_u_FQ.spectral_zeroAliasingModes();
	io_v_FQ.spectral_zeroAliasingModes();
}

#endif


