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
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_rexi");

#if !SWEET_USE_LIBFFT
	std::cerr << "Spectral space required for solvers, use compile option --libfft=enable" << std::endl;
	exit(-1);
#endif

	PDESWEPlaneTS_BaseInterface::setup(io_ops);

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
		const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
		sweet::PlaneData_Spectral &o_u,	///< prognostic variables
		sweet::PlaneData_Spectral &o_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	/// WARNING: i_h_pert might be identical to o_h_pert
	run_timestep_real(i_h_pert, i_u, i_v, o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);
}



bool SWE_Plane_Mori_Zwanzig_TS_l_exp::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	PDESWEPlaneTS_BaseInterface::shackRegistration(io_shackDict);
	ts_l_direct.shackRegistration(io_shackDict);

	return error.exists();
}

void SWE_Plane_Mori_Zwanzig_TS_l_exp::run_timestep_real(
		const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,		///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,		///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
		sweet::PlaneData_Spectral &o_u,			///< prognostic variables
		sweet::PlaneData_Spectral &o_v,			///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{

	final_timestep = false;

	if (exp_use_direct_solution)
	{

		o_h_pert = i_h_pert;
		o_u = i_u;
		o_v = i_v;
		ts_l_direct.runTimestep(o_h_pert, o_u, o_v, i_dt, i_simulation_timestamp);

		return;
	}

}


void SWE_Plane_Mori_Zwanzig_TS_l_exp::runTimestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	runTimestep(
			io_h, io_u, io_v,
			io_h, io_u, io_v,
			i_dt,
			i_simulation_timestamp
		);
}

SWE_Plane_Mori_Zwanzig_TS_l_exp::~SWE_Plane_Mori_Zwanzig_TS_l_exp()
{
}

