/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_erk.hpp"


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_erk::euler_timestep_update(
		const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	// A- grid method
	if (!shackPlaneDataOps->space_grid_use_c_staggering)
	{
		/*
		 * linearized non-conservative (advective) formulation:
		 *
		 * h_t = -h0*u_x - h0*v_ym
		 * u_t = -g * h_x + f*v
		 * v_t = -g * h_y - f*u
		 */

#if 1
		o_u_t = -shackPDESWEPlane->gravitation*ops->diff_c_x(i_h) + shackPDESWEPlane->plane_rotating_f0*i_v;
		o_v_t = -shackPDESWEPlane->gravitation*ops->diff_c_y(i_h) - shackPDESWEPlane->plane_rotating_f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_x(i_u) + ops->diff_c_y(i_v))*shackPDESWEPlane->h0;
#else

	#if 0
		// U-only
		o_u_t = -shackPDESWEPlane->gravitation*ops->diff_c_x(i_h) + shackPDESWEPlane->plane_rotating_f0*i_v;
		//o_v_t.physical_set_zero();
		o_v_t = - shackPDESWEPlane->plane_rotating_f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_x(i_u))*shackPDESWEPlane->h0;

	#else
		// V-only
		//o_u_t.spectral_set_zero();
		o_u_t = +shackPDESWEPlane->plane_rotating_f0*i_v;
		o_v_t = -shackPDESWEPlane->gravitation*ops->diff_c_y(i_h) - shackPDESWEPlane->plane_rotating_f0*i_u;// - shackDict.sim.f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_y(i_v))*shackPDESWEPlane->h0;
	#endif

#endif
	}
	else // shackDict.disc.use_staggering = true
	{
		// STAGGERED GRID

		/*
		 * Sadourny energy conserving scheme
		 *
		 * Note, that this grid does not follow the formulation
		 * in the paper of Robert Sadourny, but looks as follows:
		 *
		 *              ^
		 *              |
		 *       ______v0,1_____
		 *       |             |
		 *       |			   |
		 *       |             |
		 *  u0,0 |->  H/P0,0   |u1,0 ->
		 *(0,0.5)|			   |
		 *       |      ^      |
		 *   q0,0|______|______|
		 * (0,0)      v0,0
		 *           (0.5,0)
		 *
		 * V_t + q N x (P V) + grad( g P + 1/2 V*V) = 0
		 * P_t + div(P V) = 0
		 */

		sweet::PlaneData_Spectral H = shackPDESWEPlane->gravitation*i_h;// + 0.5*(ops->avg_f_x(i_u*i_u) + ops->avg_f_y(i_v*i_v));

		sweet::PlaneData_Physical o_u_t_phys(o_u_t.planeDataConfig);
		sweet::PlaneData_Physical o_v_t_phys(o_v_t.planeDataConfig);
		o_u_t_phys = ops->avg_f_y(shackPDESWEPlane->plane_rotating_f0*ops->avg_b_x(i_v.toPhys())) - ops->diff_b_x(H).toPhys();
		o_v_t_phys = -ops->avg_f_x(shackPDESWEPlane->plane_rotating_f0*ops->avg_b_y(i_u.toPhys())) - ops->diff_b_y(H).toPhys();
		o_u_t.loadPlaneDataPhysical(o_u_t_phys);
		o_v_t.loadPlaneDataPhysical(o_v_t_phys);

		/*
		 * P UPDATE
		 */
		o_h_t = -ops->diff_f_x(shackPDESWEPlane->h0*i_u) - ops->diff_f_y(shackPDESWEPlane->h0*i_v);
	}
}



void SWE_Plane_TS_l_erk::runTimestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_erk: Only constant time step size allowed (please set --dt)");

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&SWE_Plane_TS_l_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
bool SWE_Plane_TS_l_erk::setup(
		sweet::PlaneOperators *io_ops,
		int i_order
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	assert(i_order > 0);
	timestepping_order = i_order;

	timestepping_rk.setupBuffers(ops->planeDataConfig, timestepping_order);

	//if (shackDict.disc.use_staggering)
	//	SWEETError("Staggering not supported for l_erk");
	return true;
}

bool SWE_Plane_TS_l_erk::setup(
		sweet::PlaneOperators *io_ops
)
{
	return setup(io_ops, shackPDESWETimeDisc->timestepping_order);
}
