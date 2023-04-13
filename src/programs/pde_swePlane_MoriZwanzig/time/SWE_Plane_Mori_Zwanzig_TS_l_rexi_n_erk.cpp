/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 */

#include "SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk.hpp"

#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


bool SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{

	ts_l_rexi.shackRegistration(io_shackDict);
	PDESWEPlaneMoriZwanzigTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::euler_timestep_update_nonlinear(
		const sweet::PlaneData_Spectral &i_h_A,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_A,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_A,	///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_B,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_B,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_B,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_u_t = -i_u_A*ops->diff_c_x(i_u_B) - i_v_A*ops->diff_c_y(i_u_B);
	o_v_t = -i_u_A*ops->diff_c_x(i_v_B) - i_v_A*ops->diff_c_y(i_v_B);
	if (use_only_linear_divergence) //only nonlinear advection left to solve
		o_h_t = - (i_u_B*ops->diff_c_x(i_h_A) + i_v_B*ops->diff_c_y(i_h_A));
	else //full nonlinear equation on h
		o_h_t = -ops->diff_c_x(i_u_B*i_h_A) - ops->diff_c_y(i_v_B*i_h_A);

}


void SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::runTimestep(
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
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_rexi_n_erk: Only constant time step size allowed");

	///////////////////////////
	// solve equation for SP //
	///////////////////////////
	if (timestepping_order_nonlinear_SP == 1)
	{

		sweet::PlaneData_Spectral h_pert_N(io_h_pert_SP.planeDataConfig);
		sweet::PlaneData_Spectral u_N(io_u_SP.planeDataConfig);
		sweet::PlaneData_Spectral v_N(io_v_SP.planeDataConfig);

		ts_l_rexi.runTimestep(
				io_h_pert_SP, io_u_SP, io_v_SP,
				i_dt,
				i_simulation_timestamp
			);

		io_h_pert_SP = 1. / shackPDESWEPlane->epsilon * io_h_pert_SP;
		io_u_SP = 1. / shackPDESWEPlane->epsilon * io_u_SP;
		io_v_SP = 1. / shackPDESWEPlane->epsilon * io_v_SP;

		timestepping_rk_SP.runTimestep(
				this,
				&SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h_pert_SP, io_u_SP, io_v_SP,
				io_h_pert_SP, io_u_SP, io_v_SP,
				h_pert_N, u_N, v_N,
				i_dt,
				timestepping_order_nonlinear_SP,
				i_simulation_timestamp
			);

		this->projection->project_SP(h_pert_N, u_N, v_N);

		io_h_pert_SP += h_pert_N;
		io_u_SP += u_N;
		io_v_SP += v_N;
	}
	else if (timestepping_order_nonlinear_SP == 2)
	{
	}
	else
		SWEETError("SWE_Plane_TS_l_rexi_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");


	///////////////////////////
	// solve equation for SQ //
	///////////////////////////
	if (timestepping_order_nonlinear_SQ == 1)
	{

		sweet::PlaneData_Spectral h_pert_N_1(io_h_pert_SP.planeDataConfig);
		sweet::PlaneData_Spectral u_N_1(io_u_SP.planeDataConfig);
		sweet::PlaneData_Spectral v_N_1(io_v_SP.planeDataConfig);
		sweet::PlaneData_Spectral h_pert_N_2(io_h_pert_SP.planeDataConfig);
		sweet::PlaneData_Spectral u_N_2(io_u_SP.planeDataConfig);
		sweet::PlaneData_Spectral v_N_2(io_v_SP.planeDataConfig);

		timestepping_rk_SQ.runTimestep(
				this,
				&SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h_pert_SQ, io_u_SQ, io_v_SQ,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				h_pert_N_1, u_N_1, v_N_1,
				i_dt,
				timestepping_order_nonlinear_SQ,
				i_simulation_timestamp
			);

		this->projection->project_SQ(h_pert_N_1, u_N_1, v_N_1);

		timestepping_rk_SQ.runTimestep(
				this,
				&SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				h_pert_N_2, u_N_2, v_N_2,
				i_dt,
				timestepping_order_nonlinear_SQ,
				i_simulation_timestamp
			);

		this->projection->project_SQ(h_pert_N_2, u_N_2, v_N_2);

		io_h_pert_SQ += h_pert_N_1 + h_pert_N_2;
		io_u_SQ += u_N_1 + u_N_2;
		io_v_SQ += v_N_1 + v_N_2;
	}
	else if (timestepping_order_nonlinear_SQ == 2)
	{
	}
	else
		SWEETError("SWE_Plane_TS_l_rexi_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");


	///////////////////////////
	// solve equation for FQ //
	///////////////////////////

	if (timestepping_order_nonlinear_FQ == 1)
	{

		sweet::PlaneData_Spectral h_pert_N_1(io_h_pert_SP.planeDataConfig);
		sweet::PlaneData_Spectral u_N_1(io_u_SP.planeDataConfig);
		sweet::PlaneData_Spectral v_N_1(io_v_SP.planeDataConfig);
		sweet::PlaneData_Spectral h_pert_N_2(io_h_pert_SP.planeDataConfig);
		sweet::PlaneData_Spectral u_N_2(io_u_SP.planeDataConfig);
		sweet::PlaneData_Spectral v_N_2(io_v_SP.planeDataConfig);

		ts_l_rexi.runTimestep(
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				i_dt,
				i_simulation_timestamp
			);

		io_h_pert_FQ = 1. / shackPDESWEPlane->epsilon * io_h_pert_FQ;
		io_u_FQ = 1. / shackPDESWEPlane->epsilon * io_u_FQ;
		io_v_FQ = 1. / shackPDESWEPlane->epsilon * io_v_FQ;


		timestepping_rk_FQ.runTimestep(
				this,
				&SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h_pert_SQ, io_u_SQ, io_v_SQ,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				h_pert_N_1, u_N_1, v_N_1,
				i_dt,
				timestepping_order_nonlinear_FQ,
				i_simulation_timestamp
			);

		this->projection->project_FQ(h_pert_N_1, u_N_1, v_N_1);

		timestepping_rk_FQ.runTimestep(
				this,
				&SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				h_pert_N_2, u_N_2, v_N_2,
				i_dt,
				timestepping_order_nonlinear_FQ,
				i_simulation_timestamp
			);

		this->projection->project_FQ(h_pert_N_2, u_N_2, v_N_2);

		io_h_pert_FQ += h_pert_N_1 + h_pert_N_2;
		io_u_FQ += u_N_1 + u_N_2;
		io_v_FQ += v_N_1 + v_N_2;

	}
	else if (timestepping_order_nonlinear_FQ == 2)
	{
		//////ts_l_rexi.runTimestep(
		//////		io_h, io_u, io_v,
		//////		i_dt*0.5,
		//////		i_simulation_timestamp
		//////	);

		//////// standard time stepping
		//////timestepping_rk.runTimestep(
		//////		this,
		//////		&SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
		//////		io_h, io_u, io_v,
		//////		i_dt,
		//////		timestepping_order_nonlinear,
		//////		i_simulation_timestamp
		//////	);

		//////ts_l_rexi.runTimestep(
		//////		io_h, io_u, io_v,
		//////		i_dt*0.5,
		//////		i_simulation_timestamp
		//////	);
	}
	else
	{
		SWEETError("SWE_Plane_TS_l_rexi_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");
	}
}



/*
 * Setup
 */
bool SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk::setup(
	sweet::PlaneOperators *io_ops
)
{
	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	ts_l_rexi.setup(io_ops, "phi0");

	timestepping_order_nonlinear_SP = shackPDESWETimeDisc->timestepping_order_SP;
	timestepping_order_nonlinear_SQ = shackPDESWETimeDisc->timestepping_order_SQ;
	timestepping_order_nonlinear_FQ = shackPDESWETimeDisc->timestepping_order_FQ;

	timestepping_rk_SP.setupBuffers(ops->planeDataConfig, timestepping_order_nonlinear_SP);
	timestepping_rk_SQ.setupBuffers(ops->planeDataConfig, timestepping_order_nonlinear_SQ);
	timestepping_rk_FQ.setupBuffers(ops->planeDataConfig, timestepping_order_nonlinear_FQ);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_rexi_n_erk");

	return true;
}

