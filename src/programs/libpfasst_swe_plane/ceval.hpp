#ifndef _CEVAL_HPP_
#define _CEVAL_HPP_

#include "PlaneDataVars.hpp"
#include "PlaneDataCtx.hpp"

/*
  Right-hand-side functions called from Fortran 
*/

extern "C"
{
  // initialization of the variables (initial condition)
  void cinitial(
		PlaneDataCtx *i_ctx,
		double i_t,
		double i_dt, 
		PlaneDataVars *o_Y
		);

  // finalizes the time step when libpfasst is done 
  void cfinal(
	      PlaneDataCtx *i_ctx, 
	      PlaneDataVars *i_Y,
	      int i_nnodes,
	      int i_niters
	      );

  // computes a reference solution to check libpfasst's results
  void creference(
		  double i_t,
		  PlaneDataCtx *i_ctx,
		  PlaneDataVars *i_Y
		  );

  // evaluates the explicit piece
  void ceval_f1(
		PlaneDataVars *i_Y, 
		double i_t, 
		PlaneDataCtx *i_ctx,
		PlaneDataVars *o_F1
		);

  // evaluates the first implicit piece
  void ceval_f2 (
		 PlaneDataVars *i_Y, 
		 double i_t, 
		 PlaneDataCtx *i_ctx,
		 PlaneDataVars *o_F2 
		 );

  // solves the first implicit system
  void ccomp_f2 (
		 PlaneDataVars *io_Y, 
		 double i_t, 
		 double i_dt, 
		 PlaneDataVars *i_Rhs, 
		 PlaneDataCtx *i_ctx,
		 PlaneDataVars *o_F2 
		 );

  // evaluates the second implicit piece
  void ceval_f3 (PlaneDataVars *i_Y, 
		 double i_t, 
		 PlaneDataCtx *i_ctx,
		 PlaneDataVars *o_F3 
		 );
		 
  // solves the second implicit system
  void ccomp_f3 (PlaneDataVars *io_Y, 
		 double i_t, 
		 double i_dt, 
		 PlaneDataVars *i_Rhs,
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F3
		 );

  // applies implicit viscosity to the variables
  void capply_viscosity(PlaneDataVars *io_Y,
			double i_t, 
			double i_dt,
			int i_level, 
			PlaneDataCtx *i_ctx
			);
  
}

#endif
