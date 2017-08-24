/*
 * Burgers_Plane_TS_ln_adomian.cpp
 *
 *  Created on: 03 August 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_ln_adomian.hpp"


void Burgers_Plane_TS_ln_adomian::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{

	// setup dummy data
	PlaneData tmp(io_u.planeDataConfig);
#if SWEET_USE_PLANE_SPECTRAL_SPACE
	tmp.spectral_set_zero();
#endif
	tmp.physical_set_zero();
	/*
	 * Calculating the analytic solution to the initial condition i_prog_u
	 * with Adomian decomposition
	 */

   //TODO: implement this correctly
   if (op.diff_c_y(io_u).reduce_maxAbs()>1e-11)
      FatalError("The analyitical solution works only in 1d!");

	std::vector<PlaneData> u;
	u.resize(20,PlaneData(io_u.planeDataConfig));
	std::vector<PlaneData> A;
	A.resize(20,PlaneData(io_u.planeDataConfig));

   //Forward differencing uses dt
   double time = i_fixed_dt;

   for (int kk = 0; kk < 20; kk++) {
      A[kk] = tmp;
      u[kk] = tmp;
   }

	u[0] = io_u;

	for (int kk = 0; kk < 19; kk++) {
		for (int ii = 0; ii < kk+1; ii++) {
			for (int jj = 0; jj < ii+1; jj++) {
				if ((ii+jj)==kk)
				{
					if (ii==jj) {
						A[kk] = A[kk] + op.diff_c_x(u[ii])*u[ii];
					} else {
						A[kk] = A[kk] + op.diff_c_x(u[ii])*u[jj] + op.diff_c_x(u[jj])*u[ii];
					}
				}
			}
		}

		u[kk+1] = simVars.sim.viscosity*op.diff2_c_x(u[kk])-A[kk];
		u[kk+1] = u[kk+1] * time / (kk+1);

		io_u=io_u+u[kk+1];
      std::cout << std::endl << u[kk+1].reduce_rms() <<std::endl<<std::endl;
		if(u[kk+1].reduce_rms() < 1e-11)
			break;
	}
}



/*
 * Setup
 */
void Burgers_Plane_TS_ln_adomian::setup()
{
}


Burgers_Plane_TS_ln_adomian::Burgers_Plane_TS_ln_adomian(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup();
}



Burgers_Plane_TS_ln_adomian::~Burgers_Plane_TS_ln_adomian()
{
}

