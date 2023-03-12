/*
 * ShackDiagnostics.hpp
 *
 *  Created on: Feb 21, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKDIAGNOSTICS_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKDIAGNOSTICS_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>
#include "ShackPDESWEPlane.hpp"



/**
 * Diagnostic variables
 */
class ShackPDESWEPlaneDiagnostics	:
	public sweet::ShackInterface
{
public:
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

	void backup_reference()
	{
		ref_total_mass = total_mass;
		ref_kinetic_energy = kinetic_energy;
		ref_potential_energy = potential_energy;
		ref_total_energy = total_energy;
		ref_total_potential_enstrophy = total_potential_enstrophy;
	}

	void printProgramArguments(const std::string& i_prefix = "")
	{
		/*
		 * Doesn't exist
		 */
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		/*
		 * Doesn't exist
		 */
		return true;
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "DIAGNOSTICS:" << std::endl;
		std::cout << " + total_mass: " << total_mass << std::endl;
		std::cout << " + total_energy: " << total_energy << std::endl;
		std::cout << " + kinetic_energy: " << kinetic_energy << std::endl;
		std::cout << " + potential_energy: " << potential_energy << std::endl;
		std::cout << " + total_potential_enstrophy: " << total_potential_enstrophy << std::endl;
		std::cout << std::endl;
	}


	void update_nonstaggered_huv_to_mass_energy_enstrophy(
			sweet::PlaneOperators &ops,
			sweet::ShackPlaneDataOps *shackPlaneDataOps,
			ShackPDESWEPlane *shackPDESWEPlane,
			sweet::PlaneData_Spectral &i_prog_h, //h perturbation
			sweet::PlaneData_Spectral &i_prog_u,
			sweet::PlaneData_Spectral &i_prog_v
	)
	{
		double normalization = (shackPlaneDataOps->plane_domain_size[0]*shackPlaneDataOps->plane_domain_size[1]) /
								((double)shackPlaneDataOps->space_res_physical[0]*(double)shackPlaneDataOps->space_res_physical[1]);

		//std::cout << "Size x, sixe y" << (shackPlaneDataOps->domain_size[0]) << (shackPlaneDataOps->domain_size[1]) << std::endl;
		//std::cout << "resphysx, resphysy" << (double)shackPlaneDataOps->res_physical[0] << (double)shackPlaneDataOps->res_physical[1] << std::endl;
		//std::cout << "normal" << normalization << std::endl;

		sweet::PlaneData_Physical h_phys = i_prog_h.toPhys();
		sweet::PlaneData_Physical u_phys = i_prog_u.toPhys();
		sweet::PlaneData_Physical v_phys = i_prog_v.toPhys();

		// mass (mean depth needs to be added)
		total_mass = (h_phys+ shackPDESWEPlane->h0).physical_reduce_sum_quad() * normalization;

		// energy
		sweet::PlaneData_Physical pot_energy = (h_phys + shackPDESWEPlane->h0)*(shackPDESWEPlane->gravitation*normalization);
		sweet::PlaneData_Physical kin_energy = (h_phys + shackPDESWEPlane->h0)*(u_phys*u_phys + v_phys*v_phys)*(0.5*normalization);

		potential_energy = pot_energy.physical_reduce_sum_quad();
		kinetic_energy = kin_energy.physical_reduce_sum_quad();

		total_energy = kinetic_energy + potential_energy;

		// absolute vorticity
		sweet::PlaneData_Spectral eta = (ops.diff_c_x(i_prog_v) - ops.diff_c_y(i_prog_u) + shackPDESWEPlane->plane_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}



public:
	void update_staggered_huv_to_mass_energy_enstrophy(
			sweet::PlaneOperators &op,
			sweet::ShackPlaneDataOps *shackPlaneDataOps,
			ShackPDESWEPlane *shackPDESWEPlane,
			sweet::PlaneData_Spectral &i_prog_h,
			sweet::PlaneData_Spectral &i_prog_u,
			sweet::PlaneData_Spectral &i_prog_v
	)
	{
		double normalization = (shackPlaneDataOps->plane_domain_size[0]*shackPlaneDataOps->plane_domain_size[1]) /
								((double)shackPlaneDataOps->space_res_physical[0]*(double)shackPlaneDataOps->space_res_physical[1]);

		sweet::PlaneData_Physical h_phys = i_prog_h.toPhys();

		// mass
		total_mass = (h_phys + shackPDESWEPlane->h0).physical_reduce_sum_quad() * normalization;

		sweet::PlaneData_Physical u_phys = op.avg_b_x(i_prog_u.toPhys());
		sweet::PlaneData_Physical v_phys = op.avg_b_y(i_prog_v.toPhys());

		sweet::PlaneData_Spectral u(u_phys.planeDataConfig);
		sweet::PlaneData_Spectral v(v_phys.planeDataConfig);
		u.loadPlaneDataPhysical(u_phys);
		v.loadPlaneDataPhysical(v_phys);

		// energy
		sweet::PlaneData_Physical pot_energy = (h_phys + shackPDESWEPlane->h0)*(shackPDESWEPlane->gravitation*normalization);
		sweet::PlaneData_Physical kin_energy = (h_phys + shackPDESWEPlane->h0)*(u_phys*u_phys+v_phys*v_phys)*(0.5*normalization);

		total_energy = (pot_energy + kin_energy).physical_reduce_sum_quad();

		// total vorticity
		sweet::PlaneData_Spectral eta = (op.diff_c_x(v) - op.diff_c_y(u) + shackPDESWEPlane->plane_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toPhys().physical_reduce_sum_quad() * normalization;
	}
};




#endif
