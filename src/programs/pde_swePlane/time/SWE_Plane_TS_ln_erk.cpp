/*
 * SWE_Plane_TS_ln_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *
 *  	2017-07-13: Updated and validated by P. Peixoto.
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_ln_erk.hpp"


/*
 * Main routine for method to be used in case of finite differences
 *
 * - A-Grid with spectral or FD spatial discretizations
 * - C-Grid with energy conserving FD scheme for spatial discretizations
 *
 */
void SWE_Plane_TS_ln_erk::euler_timestep_update(
		const sweet::PlaneData_Spectral &i_h,	///< prognostic variables (perturbed part of height)
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	// A-grid method
	if (!shackPlaneDataOps->space_grid_use_c_staggering)
	{
		/*
		 * non-conservative (advective) formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
		 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
		 */

		sweet::PlaneData_Spectral total_h = i_h + shackPDESWEPlane->h0;

		o_u_t = -shackPDESWEPlane->gravitation*ops->diff_c_x(total_h) - i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
		o_v_t = -shackPDESWEPlane->gravitation*ops->diff_c_y(total_h) - i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);

		o_u_t += shackPDESWEPlane->plane_rotating_f0*i_v;
		o_v_t -= shackPDESWEPlane->plane_rotating_f0*i_u;

		// standard update
		/*
		 * P UPDATE
		 */
		if (!use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			//o_h_t = -ops->diff_f_x(U) - ops->diff_f_y(V);
			o_h_t = -ops->diff_c_x(i_u*total_h) - ops->diff_c_y(i_v*total_h);
		}
		else // use linear divergence
		{
			//o_h_t = -ops->diff_f_x(shackPDESWEPlane->h0*i_u) - ops->diff_f_y(shackPDESWEPlane->h0*i_v);
			o_h_t = -i_u*ops->diff_c_x(total_h) - i_v*ops->diff_c_y(total_h) + //nonlinear adv
					-ops->diff_c_x(i_u*shackPDESWEPlane->h0) - ops->diff_c_y(i_v*shackPDESWEPlane->h0); //linear div
		}

	}
	else // simVars.disc.use_staggering = true
	{
		// STAGGERED GRID

		sweet::PlaneData_Spectral U(i_h.planeDataConfig); // U flux
		sweet::PlaneData_Spectral V(i_h.planeDataConfig); // V flux
		sweet::PlaneData_Spectral H(i_h.planeDataConfig); //Bernoulli potential

		sweet::PlaneData_Physical U_phys(i_h.planeDataConfig); // U flux
		sweet::PlaneData_Physical V_phys(i_h.planeDataConfig); // V flux
		sweet::PlaneData_Physical H_phys(i_h.planeDataConfig); //Bernoulli potential

		sweet::PlaneData_Physical i_u_phys = i_u.toPhys();
		sweet::PlaneData_Physical i_v_phys = i_v.toPhys();

		sweet::PlaneData_Physical total_h_phys = i_h.toPhys() + shackPDESWEPlane->h0;
		sweet::PlaneData_Spectral total_h(i_h.planeDataConfig);
		total_h.loadPlaneDataPhysical(total_h_phys);


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
		/*
		 * U and V updates
		 */


		U_phys = ops->avg_b_x(total_h_phys)*i_u_phys;
		V_phys = ops->avg_b_y(total_h_phys)*i_v_phys;
		H_phys = shackPDESWEPlane->gravitation*total_h_phys + 0.5*(ops->avg_f_x(i_u_phys*i_u_phys) + ops->avg_f_y(i_v_phys*i_v_phys));

		U.loadPlaneDataPhysical(U_phys);
		V.loadPlaneDataPhysical(V_phys);
		H.loadPlaneDataPhysical(H_phys);


		// Potential vorticity
		sweet::PlaneData_Physical total_h_pv_phys = total_h_phys;
		sweet::PlaneData_Spectral total_h_pv = total_h_phys(i_h.planeDataConfig);
		total_h_pv_phys = ops->avg_b_x(ops->avg_b_y(total_h_phys));
		total_h_pv.loadPlaneDataPhysical(total_h_pv_phys);

#if 0
		if (total_h_pv.reduce_min() < 0.00000001)
		{
			std::cerr << "Test case not adequate for vector invariant formulation. Null or negative water height" << std::endl;
			std::cerr << "Min h_pv   : " << total_h_pv.reduce_min() << std::endl;
			std::cerr << "Min h_total: " << total_h.reduce_min() << std::endl;
			std::cerr << "Min h_pert : " << i_h.reduce_min() << std::endl;
			SWEETError("SWE_Plane_TS_ln_erk: Methods unstable or inadequate for vector invariant swe");;
		}
#endif

		sweet::PlaneData_Spectral q = (ops->diff_b_x(i_v) - ops->diff_b_y(i_u) + shackPDESWEPlane->plane_rotating_f0) / total_h_pv;
		sweet::PlaneData_Physical q_phys = q.toPhys();

		// u, v tendencies
		// Energy conserving scheme
		sweet::PlaneData_Physical o_u_t_phys = ops->avg_f_y(q_phys*ops->avg_b_x(V_phys)) - ops->diff_b_x(H).toPhys();
		sweet::PlaneData_Physical o_v_t_phys = -ops->avg_f_x(q_phys*ops->avg_b_y(U_phys)) - ops->diff_b_y(H).toPhys();
		o_u_t.loadPlaneDataPhysical(o_u_t_phys);
		o_v_t.loadPlaneDataPhysical(o_v_t_phys);

		/*
		 * P UPDATE
		 */
		if (!use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			o_h_t = -ops->diff_f_x(U) - ops->diff_f_y(V);
		}
		else // use linear divergence
		{
			o_h_t = -i_u*ops->diff_f_x(total_h) - i_v*ops->diff_f_y(total_h) + //nonlinear adv
					-ops->diff_f_x(i_u*shackPDESWEPlane->h0) - ops->diff_f_y(i_v*shackPDESWEPlane->h0); //linear div
			//o_h_t = -ops->diff_f_x(shackPDESWEPlane->h0*i_u) - ops->diff_f_y(shackPDESWEPlane->h0*i_v);
		}
	}
}



void SWE_Plane_TS_ln_erk::run_timestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_ln_erk: Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Plane_TS_ln_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
bool SWE_Plane_TS_ln_erk::setup(
		sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	return true;
}

