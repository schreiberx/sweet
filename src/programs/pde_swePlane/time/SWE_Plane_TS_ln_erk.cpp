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
	if (!simVars.disc.space_grid_use_c_staggering)
	{
		/*
		 * non-conservative (advective) formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
		 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
		 */

		sweet::PlaneData_Spectral total_h = i_h + simVars.sim.h0;

		o_u_t = -simVars.sim.gravitation*op.diff_c_x(total_h) - i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
		o_v_t = -simVars.sim.gravitation*op.diff_c_y(total_h) - i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

		o_u_t += simVars.sim.plane_rotating_f0*i_v;
		o_v_t -= simVars.sim.plane_rotating_f0*i_u;

		// standard update
		/*
		 * P UPDATE
		 */
		if (!use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			//o_h_t = -op.diff_f_x(U) - op.diff_f_y(V);
			o_h_t = -op.diff_c_x(i_u*total_h) - op.diff_c_y(i_v*total_h);
		}
		else // use linear divergence
		{
			//o_h_t = -op.diff_f_x(simVars.sim.h0*i_u) - op.diff_f_y(simVars.sim.h0*i_v);
			o_h_t = -i_u*op.diff_c_x(total_h) - i_v*op.diff_c_y(total_h) + //nonlinear adv
					-op.diff_c_x(i_u*simVars.sim.h0) - op.diff_c_y(i_v*simVars.sim.h0); //linear div
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

		sweet::PlaneData_Physical total_h_phys = i_h.toPhys() + simVars.sim.h0;
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


		U_phys = op.avg_b_x(total_h_phys)*i_u_phys;
		V_phys = op.avg_b_y(total_h_phys)*i_v_phys;
		H_phys = simVars.sim.gravitation*total_h_phys + 0.5*(op.avg_f_x(i_u_phys*i_u_phys) + op.avg_f_y(i_v_phys*i_v_phys));

		U.loadPlaneDataPhysical(U_phys);
		V.loadPlaneDataPhysical(V_phys);
		H.loadPlaneDataPhysical(H_phys);


		// Potential vorticity
		sweet::PlaneData_Physical total_h_pv_phys = total_h_phys;
		sweet::PlaneData_Spectral total_h_pv = total_h_phys(i_h.planeDataConfig);
		total_h_pv_phys = op.avg_b_x(op.avg_b_y(total_h_phys));
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

		sweet::PlaneData_Spectral q = (op.diff_b_x(i_v) - op.diff_b_y(i_u) + simVars.sim.plane_rotating_f0) / total_h_pv;
		sweet::PlaneData_Physical q_phys = q.toPhys();

		// u, v tendencies
		// Energy conserving scheme
		sweet::PlaneData_Physical o_u_t_phys = op.avg_f_y(q_phys*op.avg_b_x(V_phys)) - op.diff_b_x(H).toPhys();
		sweet::PlaneData_Physical o_v_t_phys = -op.avg_f_x(q_phys*op.avg_b_y(U_phys)) - op.diff_b_y(H).toPhys();
		o_u_t.loadPlaneDataPhysical(o_u_t_phys);
		o_v_t.loadPlaneDataPhysical(o_v_t_phys);

		/*
		 * P UPDATE
		 */
		if (!use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			o_h_t = -op.diff_f_x(U) - op.diff_f_y(V);
		}
		else // use linear divergence
		{
			o_h_t = -i_u*op.diff_f_x(total_h) - i_v*op.diff_f_y(total_h) + //nonlinear adv
					-op.diff_f_x(i_u*simVars.sim.h0) - op.diff_f_y(i_v*simVars.sim.h0); //linear div
			//o_h_t = -op.diff_f_x(simVars.sim.h0*i_u) - op.diff_f_y(simVars.sim.h0*i_v);
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
void SWE_Plane_TS_ln_erk::setup(
		int i_order,	///< order of RK time stepping method
		bool i_use_only_linear_divergence
)
{
	timestepping_order = i_order;
	use_only_linear_divergence = i_use_only_linear_divergence;
}


SWE_Plane_TS_ln_erk::SWE_Plane_TS_ln_erk(
		sweet::ShackDictionary *shackDict,
		sweet::PlaneOperators &i_op
)	:
		shackDict(io_shackDict),
		op(i_op)
{
	setup(simVars.disc.timestepping_order, false);
}



SWE_Plane_TS_ln_erk::~SWE_Plane_TS_ln_erk()
{
}

