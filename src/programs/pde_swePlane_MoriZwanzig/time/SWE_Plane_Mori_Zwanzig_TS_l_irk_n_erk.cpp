/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 */


#include "SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk.hpp"



bool SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	ts_l_irk.shackRegistration(io_shackDict);
	ts_n_erk.shackRegistration(io_shackDict);
	PDESWEPlaneMoriZwanzigTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}



void SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk::runTimestep(
		sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SP,		///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SP,		///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SQ,		///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SQ,		///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_FQ,		///< prognostic variables
		sweet::PlaneData_Spectral &io_v_FQ,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_irk_n_erk: Only constant time step size allowed");

	sweet::PlaneData_Spectral dummy(io_u_SP.planeDataConfig);

	//////sweet::PlaneData_Spectral h_linear_t1 = io_h;
	//////sweet::PlaneData_Spectral u_linear_t1 = io_u;
	//////sweet::PlaneData_Spectral v_linear_t1 = io_v;

	//////ts_l_irk.runTimestep(
	//////		h_linear_t1, u_linear_t1, v_linear_t1,
	//////		i_dt,
	//////		i_simulation_timestamp
	//////	);

	//////// compute non-linear tendencies at half time step
	//////sweet::PlaneData_Spectral h_dt_nonlinear(ops->planeDataConfig);
	//////sweet::PlaneData_Spectral u_dt_nonlinear(ops->planeDataConfig);
	//////sweet::PlaneData_Spectral v_dt_nonlinear(ops->planeDataConfig);

	//////// standard time stepping
	//////euler_timestep_update_nonlinear(
	//////		io_h, io_u, io_v,
	//////		h_dt_nonlinear, u_dt_nonlinear, v_dt_nonlinear
	//////	);

	//////io_h = h_linear_t1 + h_dt_nonlinear*i_dt;
	//////io_u = u_linear_t1 + u_dt_nonlinear*i_dt;
	//////io_v = v_linear_t1 + v_dt_nonlinear*i_dt;


	///////////////////////////
	// solve equation for SP //
	///////////////////////////
	if (equation == "P")
	{
		this->ts_l_irk.runTimestep(
				io_h_pert_SP, io_u_SP, io_v_SP,
				dummy, dummy, dummy,
				dummy, dummy, dummy,
				i_dt,
				i_simulation_timestamp
			);

		this->ts_n_erk.runTimestep_P(
				io_h_pert_SP, io_u_SP, io_v_SP,
				i_dt,
				i_simulation_timestamp
			);
	}
	else if (equation == "Q")
	{
		////////////////////////////////
		// solve system for SQ and FQ //
		////////////////////////////////

		this->ts_l_irk.runTimestep(
				dummy, dummy, dummy,
				dummy, dummy, dummy,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				i_dt,
				i_simulation_timestamp
			);

		this->ts_n_erk.runTimestep_Q(
				io_h_pert_SQ, io_u_SQ, io_v_SQ,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				i_dt,
				i_simulation_timestamp
			);
	}
	else if (equation == "SF")
	{
		////////////////////////////////////////////
		// solve system for S and F (full system) //
		////////////////////////////////////////////
		// Solutions : io_*_SP --> S
		//             io_*_FQ --> F

		this->ts_l_irk.runTimestep(
				io_h_pert_SP, io_u_SP, io_v_SP,
				dummy, dummy, dummy,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				i_dt,
				i_simulation_timestamp
			);

		this->ts_n_erk.runTimestep_SF(
				io_h_pert_SP, io_u_SP, io_v_SP,
				io_h_pert_FQ, io_u_FQ, io_v_FQ,
				i_dt,
				i_simulation_timestamp
			);
	}

}



/*
 * Setup
 */
bool SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk::setup(
	sweet::PlaneOperators *io_ops,
	std::string i_equation
)
{
	//PDESWEPlaneMoriZwanzigTS_BaseInterface::setup(io_ops);

	equation = i_equation;

	assert(io_ops != nullptr);
	assert(shackPDESWETimeDisc != nullptr);
	assert(shackPDESWEPlane != nullptr);

	timestepping_order_linear_P = shackPDESWETimeDisc->timestepping_order_P;
	timestepping_order_linear_Q = shackPDESWETimeDisc->timestepping_order_Q;
	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	ts_l_irk.setup(io_ops, timestepping_order_linear_P, timestepping_order_linear_Q);
	ts_n_erk.setup(io_ops, io_ops->planeDataConfig, equation);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_irk_n_erk");


	if (timestepping_order_linear_P != 1 || timestepping_order_linear_Q != 1)
		SWEETError("SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk: Only 1st order TS supported with this implementation. Please set --timestepping-order 1.");

	timestepping_order_nonlinear_P = timestepping_order_linear_P;
	timestepping_order_nonlinear_Q = timestepping_order_linear_Q;
	return true;
}

