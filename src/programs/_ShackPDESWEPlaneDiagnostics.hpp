/*
 *  Created on: 25 Feb 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_
#define SRC_INCLUDE_SWEET_PLANE_PLANEDIAGNOSTICS_HPP_

#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>



class _ShackPDESWEPlaneDiagnostics	:
		public ShackInterface
{
public:
	sweet::ShackDictionary *shackDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackPDESWEPlane *shackPDESWEPlane;


	/// total mass
	double total_mass = 0;

	/// kinetic energy
	double kinetic_energy = 0;

	/// potential energy
	double potential_energy = 0;

	/// total energy
	double total_energy = 0;

	/// total potential enstropy
	double total_potential_enstrophy = 0;

	double ref_total_mass = -1;
	double ref_kinetic_energy = -1;
	double ref_potential_energy = -1;
	double ref_total_energy = -1;
	double ref_total_potential_enstrophy = -1;


	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		return shackRegistration(io_shackDict);
	}

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackPlaneDataOps = shackDict->getAutoRegistration<ShackPlaneDataOps>();
		return true;
	}

	void update_nonstaggered_huv_to_mass_energy_enstrophy(
			PlaneOperators &ops,
			PlaneData_Spectral &i_prog_h, //h perturbation
			PlaneData_Spectral &i_prog_u,
			PlaneData_Spectral &i_prog_v
	)
	{
		double normalization = (shackPlaneDataOps->plane_domain_size[0]*shackPlaneDataOps->plane_domain_size[1]) /
								((double)shackPlaneDataOps->space_res_physical[0]*(double)shackPlaneDataOps->space_res_physical[1]);

		//std::cout << "Size x, sixe y" << (shackPlaneDataOps->domain_size[0]) << (shackPlaneDataOps->domain_size[1]) << std::endl;
		//std::cout << "resphysx, resphysy" << (double)shackPlaneDataOps->res_physical[0] << (double)shackPlaneDataOps->res_physical[1] << std::endl;
		//std::cout << "normal" << normalization << std::endl;

		PlaneData_Physical h_phys = i_prog_h.toPhys();
		PlaneData_Physical u_phys = i_prog_u.toPhys();
		PlaneData_Physical v_phys = i_prog_v.toPhys();

		// mass (mean depth needs to be added)
		total_mass = (h_phys+ shackPDESWEPlane->h0).physical_reduce_sum_quad() * normalization;

		// energy
		PlaneData_Physical pot_energy = (h_phys + shackPlaneDataOps->h0)*(shackPlaneDataOps->gravitation*normalization);
		PlaneData_Physical kin_energy = (h_phys + shackPlaneDataOps->h0)*(u_phys*u_phys + v_phys*v_phys)*(0.5*normalization);

		potential_energy = pot_energy.physical_reduce_sum_quad();
		kinetic_energy = kin_energy.physical_reduce_sum_quad();

		total_energy = kinetic_energy + potential_energy;

		// absolute vorticity
		PlaneData_Spectral eta = (ops.diff_c_x(i_prog_v) - ops.diff_c_y(i_prog_u) + shackPlaneDataOps->plane_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}



public:
	void update_staggered_huv_to_mass_energy_enstrophy(
			PlaneOperators &op,
			PlaneData_Spectral &i_prog_h,
			PlaneData_Spectral &i_prog_u,
			PlaneData_Spectral &i_prog_v
	)
	{
		double normalization = (shackPlaneDataOps->plane_domain_size[0]*shackPlaneDataOps->plane_domain_size[1]) /
								((double)shackPlaneDataOps->space_res_physical[0]*(double)shackPlaneDataOps->space_res_physical[1]);

		PlaneData_Physical h_phys = i_prog_h.toPhys();

		// mass
		total_mass = (h_phys + shackPlaneDataOps->h0).physical_reduce_sum_quad() * normalization;

		PlaneData_Physical u_phys = op.avg_b_x(i_prog_u.toPhys());
		PlaneData_Physical v_phys = op.avg_b_y(i_prog_v.toPhys());

		PlaneData_Spectral u(u_phys.planeDataConfig);
		PlaneData_Spectral v(v_phys.planeDataConfig);
		u.loadPlaneDataPhysical(u_phys);
		v.loadPlaneDataPhysical(v_phys);

		// energy
		PlaneData_Physical pot_energy = (h_phys + shackPlaneDataOps->h0)*(shackPlaneDataOps->gravitation*normalization);
		PlaneData_Physical kin_energy = (h_phys + shackPlaneDataOps->h0)*(u_phys*u_phys+v_phys*v_phys)*(0.5*normalization);

		total_energy = (pot_energy + kin_energy).physical_reduce_sum_quad();

		// total vorticity
		PlaneData_Spectral eta = (op.diff_c_x(v) - op.diff_c_y(u) + shackPlaneDataOps->plane_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}


	void update_nonstaggered_huv_to_mass_energy_enstrophy(
			PlaneOperators &ops,
			PlaneData_Spectral &i_prog_h, //h perturbation
			PlaneData_Spectral &i_prog_u,
			PlaneData_Spectral &i_prog_v
	)
	{
		double normalization = (shackPlaneDataOps->plane_domain_size[0]*shackPlaneDataOps->plane_domain_size[1]) /
								((double)shackPlaneDataOps->space_res_physical[0]*(double)shackPlaneDataOps->space_res_physical[1]);

		//std::cout << "Size x, sixe y" << (shackPlaneDataOps->domain_size[0]) << (shackPlaneDataOps->domain_size[1]) << std::endl;
		//std::cout << "resphysx, resphysy" << (double)shackPlaneDataOps->res_physical[0] << (double)shackPlaneDataOps->res_physical[1] << std::endl;
		//std::cout << "normal" << normalization << std::endl;

		PlaneData_Physical h_phys = i_prog_h.toPhys();
		PlaneData_Physical u_phys = i_prog_u.toPhys();
		PlaneData_Physical v_phys = i_prog_v.toPhys();

		// mass (mean depth needs to be added)
		total_mass = (h_phys+ shackPDESWEPlane->h0).physical_reduce_sum_quad() * normalization;

		// energy
		PlaneData_Physical pot_energy = (h_phys + shackPlaneDataOps->h0)*(shackPlaneDataOps->gravitation*normalization);
		PlaneData_Physical kin_energy = (h_phys + shackPlaneDataOps->h0)*(u_phys*u_phys + v_phys*v_phys)*(0.5*normalization);

		potential_energy = pot_energy.physical_reduce_sum_quad();
		kinetic_energy = kin_energy.physical_reduce_sum_quad();

		total_energy = kinetic_energy + potential_energy;

		// absolute vorticity
		PlaneData_Spectral eta = (ops.diff_c_x(i_prog_v) - ops.diff_c_y(i_prog_u) + shackPlaneDataOps->plane_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}



public:
	void update_staggered_huv_to_mass_energy_enstrophy(
			PlaneOperators &op,
			PlaneData_Spectral &i_prog_h,
			PlaneData_Spectral &i_prog_u,
			PlaneData_Spectral &i_prog_v
	)
	{
		double normalization = (shackPlaneDataOps->plane_domain_size[0]*shackPlaneDataOps->plane_domain_size[1]) /
								((double)shackPlaneDataOps->space_res_physical[0]*(double)shackPlaneDataOps->space_res_physical[1]);

		PlaneData_Physical h_phys = i_prog_h.toPhys();

		// mass
		total_mass = (h_phys + shackPlaneDataOps->h0).physical_reduce_sum_quad() * normalization;

		PlaneData_Physical u_phys = op.avg_b_x(i_prog_u.toPhys());
		PlaneData_Physical v_phys = op.avg_b_y(i_prog_v.toPhys());

		PlaneData_Spectral u(u_phys.planeDataConfig);
		PlaneData_Spectral v(v_phys.planeDataConfig);
		u.loadPlaneDataPhysical(u_phys);
		v.loadPlaneDataPhysical(v_phys);

		// energy
		PlaneData_Physical pot_energy = (h_phys + shackPlaneDataOps->h0)*(shackPlaneDataOps->gravitation*normalization);
		PlaneData_Physical kin_energy = (h_phys + shackPlaneDataOps->h0)*(u_phys*u_phys+v_phys*v_phys)*(0.5*normalization);

		total_energy = (pot_energy + kin_energy).physical_reduce_sum_quad();

		// total vorticity
		PlaneData_Spectral eta = (op.diff_c_x(v) - op.diff_c_y(u) + shackPlaneDataOps->plane_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}
};

#endif
