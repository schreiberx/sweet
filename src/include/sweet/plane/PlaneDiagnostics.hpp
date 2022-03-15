/*
 * PlaneDiagnostics.hpp
 *
 *  Created on: 25 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_
#define SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_



class PlaneDiagnostics
{
public:
	static

	void update_nonstaggered_huv_to_mass_energy_enstrophy(
			PlaneOperators &op,
			PlaneData_Spectral &i_prog_h, //h perturbation
			PlaneData_Spectral &i_prog_u,
			PlaneData_Spectral &i_prog_v,
			SimulationVariables &io_simVars
	)
	{
		double normalization = (io_simVars.sim.plane_domain_size[0]*io_simVars.sim.plane_domain_size[1]) /
								((double)io_simVars.disc.space_res_physical[0]*(double)io_simVars.disc.space_res_physical[1]);

		//std::cout << "Size x, sixe y" << (io_simVars.sim.domain_size[0]) << (io_simVars.sim.domain_size[1]) << std::endl;
		//std::cout << "resphysx, resphysy" << (double)io_simVars.disc.res_physical[0] << (double)io_simVars.disc.res_physical[1] << std::endl;
		//std::cout << "normal" << normalization << std::endl;

		PlaneData_Physical h_phys = i_prog_h.toPhys();
		PlaneData_Physical u_phys = i_prog_u.toPhys();
		PlaneData_Physical v_phys = i_prog_v.toPhys();

		// mass (mean depth needs to be added)
		io_simVars.diag.total_mass = (i_prog_h+ io_simVars.sim.h0).reduce_sum_quad() * normalization;

		// energy
		PlaneData_Physical pot_energy = (h_phys + io_simVars.sim.h0)*(io_simVars.sim.gravitation*normalization);
		PlaneData_Physical kin_energy = (h_phys + io_simVars.sim.h0)*(u_phys*u_phys + v_phys*v_phys)*(0.5*normalization);

		io_simVars.diag.potential_energy = pot_energy.physical_reduce_sum_quad();
		io_simVars.diag.kinetic_energy = kin_energy.physical_reduce_sum_quad();

		io_simVars.diag.total_energy = io_simVars.diag.kinetic_energy + io_simVars.diag.potential_energy;

		// absolute vorticity
		PlaneData_Spectral eta = (op.diff_c_x(i_prog_v) - op.diff_c_y(i_prog_u) + io_simVars.sim.plane_rotating_f0);

		// enstrophy
		io_simVars.diag.total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}



public:
	static
	void update_staggered_huv_to_mass_energy_enstrophy(
			PlaneOperators &op,
			PlaneData_Spectral &i_prog_h,
			PlaneData_Spectral &i_prog_u,
			PlaneData_Spectral &i_prog_v,
			SimulationVariables &io_simVars
	)
	{
		double normalization = (io_simVars.sim.plane_domain_size[0]*io_simVars.sim.plane_domain_size[1]) /
								((double)io_simVars.disc.space_res_physical[0]*(double)io_simVars.disc.space_res_physical[1]);

		PlaneData_Physical h_phys = i_prog_h.toPhys();

		// mass
		io_simVars.diag.total_mass = (h_phys + io_simVars.sim.h0).physical_reduce_sum_quad() * normalization;

		PlaneData_Spectral u = op.avg_b_x(i_prog_u);
		PlaneData_Spectral v = op.avg_b_y(i_prog_v);

		PlaneData_Physical u_phys = u.toPhys();
		PlaneData_Physical v_phys = v.toPhys();

		// energy
		PlaneData_Physical pot_energy = (h_phys + io_simVars.sim.h0)*(io_simVars.sim.gravitation*normalization);
		PlaneData_Physical kin_energy = (h_phys + io_simVars.sim.h0)*(u_phys*u_phys+v_phys*v_phys)*(0.5*normalization);

		io_simVars.diag.total_energy = (pot_energy + kin_energy).toPhys().physical_reduce_sum_quad();

		// total vorticity
		PlaneData_Spectral eta = (op.diff_c_x(v) - op.diff_c_y(u) + io_simVars.sim.plane_rotating_f0);

		// enstrophy
		io_simVars.diag.total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_ */
