#ifndef _CEVAL_HPP_
#define _CEVAL_HPP_

#include "SphereDataVars.hpp"
#include "SphereDataCtx.hpp"

/**
 * Write file to data and return string of file name
 */

std::string write_file(
		       SphereDataCtx  &i_ctx,
		       const SphereDataSpectral &i_sphereData,
		       const char* i_name	///< name of output variable
		       );

/*
  Right-hand-side functions called from Fortran 
*/

extern "C"
{
  // initialization of the variables (initial condition)
  void cinitial(
		SphereDataCtx *i_ctx,
		double i_t,
		double i_dt, 
		SphereDataVars *o_Y
		);

  // finalizes the time step when libpfasst is done 
  void cfinal(
	      SphereDataCtx *i_ctx, 
	      SphereDataVars *i_Y,
	      int i_nnodes,
	      int i_niters
	      );

  // evaluates the explicit piece
  void ceval_f1(
		SphereDataVars *i_Y, 
		double i_t, 
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F1
		);

  // evaluates the first implicit piece
  void ceval_f2 (
		 SphereDataVars *i_Y, 
		 double i_t, 
		 SphereDataCtx *i_ctx,
		 SphereDataVars *o_F2 
		 );

  // solves the first implicit system
  void ccomp_f2 (
		 SphereDataVars *io_Y, 
		 double i_t, 
		 double i_dt, 
		 SphereDataVars *i_Rhs, 
		 SphereDataCtx *i_ctx,
		 SphereDataVars *o_F2 
		 );

  // evaluates the second implicit piece
  void ceval_f3 (SphereDataVars *i_Y, 
		 double i_t, 
		 int i_level,
		 SphereDataCtx *i_ctx,
		 SphereDataVars *o_F3 
		 );
		 
  // solves the second implicit system
  void ccomp_f3 (SphereDataVars *io_Y, 
		 double i_t, 
		 double i_dt,
		 int i_level,
		 SphereDataVars *i_Rhs,
		 SphereDataCtx *i_ctx, 
		 SphereDataVars *o_F3
		 );

  // applies artificial diffusion
  void cfinalize (SphereDataVars *io_Y,
		  double i_t,
		  double i_dt,
		  SphereDataCtx *i_ctx);
  
}

#endif
