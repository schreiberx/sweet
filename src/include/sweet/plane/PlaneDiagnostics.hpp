/*
 * PlaneDiagnostics.hpp
 *
 *  Created on: 25 Feb 2017
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_
#define SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_



class PlaneDiagnostics
{
public:
	static
	void update_nonstaggered_h_u_v(
			PlaneOperators &op,
			PlaneData &i_prog_h,
			PlaneData &i_prog_u,
			PlaneData &i_prog_v,
			SimulationVariables &io_simVars
	)
	{
		double normalization = (io_simVars.sim.domain_size[0]*io_simVars.sim.domain_size[1]) /
								((double)io_simVars.disc.res_physical[0]*(double)io_simVars.disc.res_physical[1]);

		// mass
		io_simVars.diag.total_mass = i_prog_h.reduce_sum_quad() * normalization;

		// energy
		PlaneData pot_energy = i_prog_h*(io_simVars.sim.gravitation*normalization);
		PlaneData kin_energy = i_prog_h*(i_prog_u*i_prog_u+i_prog_v*i_prog_v)*(0.5*normalization);

		io_simVars.diag.total_energy = (pot_energy + kin_energy).reduce_sum_quad();

		// total vorticity
		PlaneData eta = (op.diff_c_x(i_prog_v) - op.diff_c_y(i_prog_u) + io_simVars.sim.f0);

		// enstrophy
		io_simVars.diag.total_potential_enstrophy = 0.5*(eta*eta).reduce_sum_quad() * normalization;
	}



public:
	static
	void update_staggered_h_u_v(
			PlaneOperators &op,
			PlaneData &i_prog_h,
			PlaneData &i_prog_u,
			PlaneData &i_prog_v,
			SimulationVariables &io_simVars
	)
	{
		double normalization = (io_simVars.sim.domain_size[0]*io_simVars.sim.domain_size[1]) /
								((double)io_simVars.disc.res_physical[0]*(double)io_simVars.disc.res_physical[1]);

		// mass
		io_simVars.diag.total_mass = i_prog_h.reduce_sum_quad() * normalization;

		PlaneData u = op.avg_b_x(i_prog_u);
		PlaneData v = op.avg_b_y(i_prog_v);

		// energy
		PlaneData pot_energy = i_prog_h*(io_simVars.sim.gravitation*normalization);
		PlaneData kin_energy = i_prog_h*(u*u+v*v)*(0.5*normalization);

		io_simVars.diag.total_energy = (pot_energy + kin_energy).reduce_sum_quad();

		// total vorticity
		PlaneData eta = (op.diff_c_x(v) - op.diff_c_y(u) + io_simVars.sim.f0);

		// enstrophy
		io_simVars.diag.total_potential_enstrophy = 0.5*(eta*eta).reduce_sum_quad() * normalization;
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_ */
