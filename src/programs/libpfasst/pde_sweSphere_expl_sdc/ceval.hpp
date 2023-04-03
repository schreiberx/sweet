#ifndef _CEVAL_HPP_
#define _CEVAL_HPP_

#include "../interface/SphereDataVars.hpp"
#include "SphereDataCtxSDC.hpp"

/**
 * Determine if output should be written for current time step & iteration
 */
bool timestep_check_output(SphereDataCtxSDC *i_ctx,
						   int i_current_iter,
						   int i_niters);

/**
 * Write file to data and return string of file name
 */

std::string write_file(
	SphereDataCtxSDC *i_ctx,
	const sweet::SphereData_Spectral &i_sphereData,
	const char *i_name ///< name of output variable
);

/*
  Right-hand-side functions called from Fortran
*/

extern "C"
{
	// initialization of the variables (initial condition)
	void cinitial(
		SphereDataCtxSDC *i_ctx,
		double i_t,
		double i_dt,
		SphereDataVars *o_Y);

	// finalizes the time step when libpfasst is done
	void cfinal(
		SphereDataCtxSDC *i_ctx,
		SphereDataVars *i_Y,
		int i_nnodes,
		int i_niters);

	// evaluates the explicit piece
	void ceval(
		SphereDataVars *i_Y,
		double i_t,
		SphereDataCtxSDC *i_ctx,
		SphereDataVars *o_F1);
}

#endif
