/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 */

#include "SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk.hpp"

#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


bool SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{

	ts_l_exp.shackRegistration(io_shackDict);
	ts_n_erk.shackRegistration(io_shackDict);
	PDESWEPlaneMoriZwanzigTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}

void SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk::runTimestep(
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
		SWEETError("SWE_Plane_TS_l_exp_n_erk: Only constant time step size allowed");

	sweet::PlaneData_Spectral dummy(io_u_SP.planeDataConfig);

	///////////////////////////
	// solve equation for SP //
	///////////////////////////
	if (equation == "P")
	{
		if (timestepping_order_nonlinear_P == 1)
		{

			ts_l_exp.runTimestep(
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
		else if (timestepping_order_nonlinear_P == 2)
		{

			ts_l_exp.runTimestep(
					io_h_pert_SP, io_u_SP, io_v_SP,
					dummy, dummy, dummy,
					dummy, dummy, dummy,
					i_dt * .5,
					i_simulation_timestamp
				);

			this->ts_n_erk.runTimestep_P(
					io_h_pert_SP, io_u_SP, io_v_SP,
					i_dt,
					i_simulation_timestamp
				);

			ts_l_exp.runTimestep(
					io_h_pert_SP, io_u_SP, io_v_SP,
					dummy, dummy, dummy,
					dummy, dummy, dummy,
					i_dt * .5,
					i_simulation_timestamp
				);

		}
		else
			SWEETError("SWE_Plane_TS_l_exp_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");
	}
	else if (equation == "Q")
	{
		////////////////////////////////
		// solve system for SQ and FQ //
		////////////////////////////////
		if (timestepping_order_nonlinear_Q == 1)
		{

			ts_l_exp.runTimestep(
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
		else if (timestepping_order_nonlinear_Q == 2)
		{

			ts_l_exp.runTimestep(
					dummy, dummy, dummy,
					dummy, dummy, dummy,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt * .5,
					i_simulation_timestamp
				);

			this->ts_n_erk.runTimestep_Q(
					io_h_pert_SQ, io_u_SQ, io_v_SQ,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt,
					i_simulation_timestamp
				);

			ts_l_exp.runTimestep(
					dummy, dummy, dummy,
					dummy, dummy, dummy,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt * .5,
					i_simulation_timestamp
				);

		}
		else
			SWEETError("SWE_Plane_TS_l_exp_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");
	}
	else if (equation == "SF")
	{
		////////////////////////////////////////////
		// solve system for S and F (full system) //
		////////////////////////////////////////////
		// Solutions : io_*_SP --> S
		//             io_*_FQ --> F
		if (timestepping_order_nonlinear_Q == 1)
		{

			ts_l_exp.runTimestep(
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
		else if (timestepping_order_nonlinear_Q == 2)
		{

			ts_l_exp.runTimestep(
					io_h_pert_SP, io_u_SP, io_v_SP,
					dummy, dummy, dummy,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt * .5,
					i_simulation_timestamp
				);

			this->ts_n_erk.runTimestep_SF(
					io_h_pert_SP, io_u_SP, io_v_SP,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt,
					i_simulation_timestamp
				);

			ts_l_exp.runTimestep(
					io_h_pert_SP, io_u_SP, io_v_SP,
					dummy, dummy, dummy,
					io_h_pert_FQ, io_u_FQ, io_v_FQ,
					i_dt * .5,
					i_simulation_timestamp
				);

		}
		else
			SWEETError("SWE_Plane_TS_l_exp_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");

	}

}



/*
 * Setup
 */
bool SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk::setup(
	sweet::PlaneOperators *io_ops,
	std::string i_equation
)
{
	ops = io_ops;
	equation = i_equation;

	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	ts_l_exp.setup(io_ops, "phi0");
	ts_n_erk.setup(io_ops, io_ops->planeDataConfig, equation);

	timestepping_order_nonlinear_P = shackPDESWETimeDisc->timestepping_order_P;
	timestepping_order_nonlinear_Q = shackPDESWETimeDisc->timestepping_order_Q;

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_exp_n_erk");

	return true;
}

