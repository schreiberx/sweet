/*
 * Burgers_Plane_TS_n_adomian.cpp
 *
 *  Created on: 03 August 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_n_adomian.hpp"


void Burgers_Plane_TS_n_adomian::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Burgers_Plane_TS_n_adomian: Only constant time step size allowed");

	// setup dummy data
	PlaneData tmp(io_u.planeDataConfig);
#if SWEET_USE_PLANE_SPECTRAL_SPACE
	tmp.spectral_set_zero();
#endif
	tmp.physical_set_zero();
	/*
	 * Applying the Adomian decomposition method. Using up to timestepping_order Adomian polynomials
	 */

	//TODO: implement this correctly
	if (op.diff_c_y(io_u).reduce_maxAbs()>1e-11)
		FatalError("The Adomian decomposition method is implemented in 1D only!");

	std::vector<PlaneData> u;
	u.resize(timestepping_order+1,PlaneData(io_u.planeDataConfig));
	std::vector<PlaneData> A;
	A.resize(timestepping_order+1,PlaneData(io_u.planeDataConfig));

	//Forward differencing uses dt
	double time = i_fixed_dt;

	for (int kk = 0; kk < timestepping_order+1; kk++) {
		A[kk] = tmp;
		u[kk] = tmp;
	}

	u[0] = io_u;

	for (int kk = 0; kk < timestepping_order; kk++) {
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

		u[kk+1] = -A[kk];
		// Integration by time, works like this, since u0 is independent of t.
		// Therefore, each integration goes from 1 to 1*t, to t*t/2, to t^2*t/3 etc.
		u[kk+1] = u[kk+1] * time / (kk+1);

		io_u=io_u+u[kk+1];

		if(u[kk+1].reduce_rms() < 1e-12)
			break;
	}
}



/*
 * Setup
 */
void Burgers_Plane_TS_n_adomian::setup(
		int i_order	///< number used of Adomian polynomials
)
{
	timestepping_order = i_order;
}


Burgers_Plane_TS_n_adomian::Burgers_Plane_TS_n_adomian(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Burgers_Plane_TS_n_adomian::~Burgers_Plane_TS_n_adomian()
{
}

