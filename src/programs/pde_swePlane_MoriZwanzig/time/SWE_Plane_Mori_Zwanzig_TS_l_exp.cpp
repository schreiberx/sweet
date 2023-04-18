/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "SWE_Plane_Mori_Zwanzig_TS_l_exp.hpp"

#include <sweet/expIntegration/REXI.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneOperatorsComplex.hpp>

#include <sweet/core/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>

#include <sweet/core/StopwatchBox.hpp>

#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
	#include <omp.h>
#endif


bool SWE_Plane_Mori_Zwanzig_TS_l_exp::setup(
	sweet::PlaneOperators *io_ops
)
{
	return setup(io_ops, "phi0");
}


bool SWE_Plane_Mori_Zwanzig_TS_l_exp::setup(
	sweet::PlaneOperators *io_ops,
	const std::string &i_function_name
)
{
	ops = io_ops;

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_rexi");

#if !SWEET_USE_LIBFFT
	std::cerr << "Spectral space required for solvers, use compile option --libfft=enable" << std::endl;
	exit(-1);
#endif

	assert(shackPlaneDataOps != nullptr);
	assert(shackExpIntegration != nullptr);

	exp_use_direct_solution = (shackExpIntegration->exp_method == "direct");

	if (exp_use_direct_solution)
	{
		ts_l_direct.setup(io_ops, i_function_name);
		return true;
	}
	else
		SWEETError("only direct exponential integration is available!");


	return true;
}


void SWE_Plane_Mori_Zwanzig_TS_l_exp::runTimestep(
/////		const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
/////		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
/////		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables
/////
/////		sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
/////		sweet::PlaneData_Spectral &o_u,	///< prognostic variables
/////		sweet::PlaneData_Spectral &o_v,	///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_pert_SP,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_SP,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_SP,		///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_pert_SQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_SQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_SQ,		///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_pert_FQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_FQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_FQ,		///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_SP,		///< prognostic variables
		sweet::PlaneData_Spectral &o_u_SP,		///< prognostic variables
		sweet::PlaneData_Spectral &o_v_SP,		///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_SQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_u_SQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_v_SQ,		///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_FQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_u_FQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_v_FQ,		///< prognostic variables


		double i_dt,
		double i_simulation_timestamp
)
{
	/// WARNING: i_h_pert might be identical to o_h_pert
	run_timestep_real(
				i_h_pert_SP, i_u_SP, i_v_SP,
				i_h_pert_SQ, i_u_SQ, i_v_SQ,
				i_h_pert_FQ, i_u_FQ, i_v_FQ,
				o_h_pert_SP, o_u_SP, o_v_SP,
				o_h_pert_SQ, o_u_SQ, o_v_SQ,
				o_h_pert_FQ, o_u_FQ, o_v_FQ,
				i_dt, i_simulation_timestamp
		);
}



bool SWE_Plane_Mori_Zwanzig_TS_l_exp::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	PDESWEPlaneMoriZwanzigTS_BaseInterface::shackRegistration(io_shackDict);
	ts_l_direct.shackRegistration(io_shackDict);

	return error.exists();
}

void SWE_Plane_Mori_Zwanzig_TS_l_exp::run_timestep_real(
//		const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
//		const sweet::PlaneData_Spectral &i_u,		///< prognostic variables
//		const sweet::PlaneData_Spectral &i_v,		///< prognostic variables
//
//		sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
//		sweet::PlaneData_Spectral &o_u,			///< prognostic variables
//		sweet::PlaneData_Spectral &o_v,			///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_pert_SP,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_SP,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_SP,		///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_pert_SQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_SQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_SQ,		///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_pert_FQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_FQ,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_FQ,		///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_SP,		///< prognostic variables
		sweet::PlaneData_Spectral &o_u_SP,		///< prognostic variables
		sweet::PlaneData_Spectral &o_v_SP,		///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_SQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_u_SQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_v_SQ,		///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_FQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_u_FQ,		///< prognostic variables
		sweet::PlaneData_Spectral &o_v_FQ,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{

	final_timestep = false;

	if (exp_use_direct_solution)
	{

		o_h_pert_SP = i_h_pert_SP;
		o_u_SP = i_u_SP;
		o_v_SP = i_v_SP;
		o_h_pert_SQ = i_h_pert_SQ;
		o_u_SQ = i_u_SQ;
		o_v_SQ = i_v_SQ;
		o_h_pert_FQ = i_h_pert_FQ;
		o_u_FQ = i_u_FQ;
		o_v_FQ = i_v_FQ;

		ts_l_direct.runTimestep(
				o_h_pert_SP, o_u_SP, o_v_SP,
				o_h_pert_SQ, o_u_SQ, o_v_SQ,
				o_h_pert_FQ, o_u_FQ, o_v_FQ,
				i_dt, i_simulation_timestamp
			);

		return;
	}

}


void SWE_Plane_Mori_Zwanzig_TS_l_exp::runTimestep(
///		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
///		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
///		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

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
	runTimestep(
			io_h_pert_SP, io_u_SP, io_v_SP,
			io_h_pert_SQ, io_u_SQ, io_v_SQ,
			io_h_pert_FQ, io_u_FQ, io_v_FQ,
			io_h_pert_SP, io_u_SP, io_v_SP,
			io_h_pert_SQ, io_u_SQ, io_v_SQ,
			io_h_pert_FQ, io_u_FQ, io_v_FQ,
			i_dt,
			i_simulation_timestamp
		);
}

SWE_Plane_Mori_Zwanzig_TS_l_exp::~SWE_Plane_Mori_Zwanzig_TS_l_exp()
{
}

