/*
 * SphereDiagnostics.hpp
 *
 *  Created on: 25 Feb 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_SPHEREDIAGNOSTICS_HPP_
#define SRC_INCLUDE_SWEET_PLANE_SPHEREDIAGNOSTICS_HPP_



class SphereDiagnostics
{
public:
	static
	void update_nonstaggered_huv_to_mass_energy_enstrophy(
			const SphereOperators &op,
			const SphereData &i_prog_h,
			const SphereData &i_prog_u,
			const SphereData &i_prog_v,
			SimulationVariables &io_simVars
	)
	{
		double normalization = (io_simVars.sim.domain_size[0]*io_simVars.sim.domain_size[1]) /
								((double)io_simVars.disc.res_physical[0]*(double)io_simVars.disc.res_physical[1]);

		// mass
		io_simVars.diag.total_mass = i_prog_h.reduce_sum_quad() * normalization;

		// energy
		SphereData pot_energy = i_prog_h*(io_simVars.sim.gravitation*normalization);
		SphereData kin_energy = i_prog_h*(i_prog_u*i_prog_u+i_prog_v*i_prog_v)*(0.5*normalization);

		io_simVars.diag.total_energy = (pot_energy + kin_energy).reduce_sum_quad();

		// total vorticity
		SphereData eta = (op.diff_c_x(i_prog_v) - op.diff_c_y(i_prog_u) + io_simVars.sim.f0);

		// enstrophy
		io_simVars.diag.total_potential_enstrophy = 0.5*(eta*eta).reduce_sum_quad() * normalization;
	}



public:
	static
	void update_staggered_huv_to_mass_energy_enstrophy(
			PlaneOperators &op,
			SphereData &i_prog_h,
			SphereData &i_prog_u,
			SphereData &i_prog_v,
			SimulationVariables &io_simVars
	)
	{
		double normalization = (io_simVars.sim.domain_size[0]*io_simVars.sim.domain_size[1]) /
								((double)io_simVars.disc.res_physical[0]*(double)io_simVars.disc.res_physical[1]);

		// mass
		io_simVars.diag.total_mass = i_prog_h.reduce_sum_quad() * normalization;

		SphereData u = op.avg_b_x(i_prog_u);
		SphereData v = op.avg_b_y(i_prog_v);

		// energy
		SphereData pot_energy = i_prog_h*(io_simVars.sim.gravitation*normalization);
		SphereData kin_energy = i_prog_h*(u*u+v*v)*(0.5*normalization);

		io_simVars.diag.total_energy = (pot_energy + kin_energy).reduce_sum_quad();

		// total vorticity
		SphereData eta = (op.diff_c_x(v) - op.diff_c_y(u) + io_simVars.sim.f0);

		// enstrophy
		io_simVars.diag.total_potential_enstrophy = 0.5*(eta*eta).reduce_sum_quad() * normalization;
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_SPHEREDIAGNOSTICS_HPP_ */
