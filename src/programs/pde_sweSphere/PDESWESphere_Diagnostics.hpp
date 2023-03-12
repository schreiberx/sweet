/*
 * PDESWESphereDiagnostics.hpp
 *
 *  Created on: Mar 10, 2023
 * Author: martin
 */

#ifndef SRC_PROGRAMS_PDE_SWESPHERE_PDESWESPHEREDIAGNOSTICS_HPP_
#define SRC_PROGRAMS_PDE_SWESPHERE_PDESWESPHEREDIAGNOSTICS_HPP_


#include <sweet/core/sphere/SphereHelpers_Integral.hpp>
#include "ShackPDESWESphere.hpp"

class PDESWESphere_Diagnostics
{
	sweet::SphereOperators *sphereOperators;
	ShackPDESWESphere *shackPDESWESphere;

	sweet::SphereHelpers_Integral sphereHelpers_integral;
	sweet::SphereData_Physical fg;

public:
	double total_mass;
	double potential_energy;
	double kinetic_energy;
	double total_potential_enstrophy;
	double total_energy;

	double ref_total_mass;
	double ref_potential_energy;
	double ref_kinetic_energy;

	int last_update_timestep_nr;

public:
	PDESWESphere_Diagnostics()	:
		sphereOperators(nullptr),
		shackPDESWESphere(nullptr),
		total_mass(0),
		potential_energy(0),
		kinetic_energy(0),
		total_potential_enstrophy(0),
		total_energy(0),

		ref_total_mass(0),
		ref_potential_energy(0),
		ref_kinetic_energy(0)
	{
	}


	void setup(
			sweet::SphereOperators *i_sphereOperators,
			ShackPDESWESphere *i_shackPDESWESphere,
			int i_verbose = 1
	)
	{
		sphereOperators = i_sphereOperators;
		shackPDESWESphere = i_shackPDESWESphere;

		sphereHelpers_integral.setup(sphereOperators->sphereDataConfig, i_verbose);

		setupFG();

		last_update_timestep_nr = -1;
	}


	bool setupFG()
	{
		if (shackPDESWESphere->sphere_use_fsphere)
			fg = sphereOperators->getFG_fSphere(shackPDESWESphere->sphere_fsphere_f0);
		else
			fg = sphereOperators->getFG_rotatingSphere(shackPDESWESphere->sphere_rotating_coriolis_omega);

		return true;
	}


public:
	void update_phi_vrt_div_2_mass_energy_enstrophy(
			const sweet::SphereOperators *i_ops,
			const sweet::SphereData_Spectral &i_prog_phi,
			const sweet::SphereData_Spectral &i_prog_vort,
			const sweet::SphereData_Spectral &i_prog_div,

			double i_sphere_radius,
			double i_gravitation
	)
	{
		assert(sphereOperators != nullptr);

		sweet::SphereData_Physical h(sphereOperators->sphereDataConfig);
		sweet::SphereData_Physical u(sphereOperators->sphereDataConfig);
		sweet::SphereData_Physical v(sphereOperators->sphereDataConfig);

		h = i_prog_phi.toPhys()*(1.0/i_gravitation);
		i_ops->vrtdiv_to_uv(i_prog_vort, i_prog_div, u, v);

		double normalization = i_sphere_radius*i_sphere_radius;

		// mass
		total_mass = sphereHelpers_integral.compute_zylinder_integral(h) * normalization;

		// energy
		//SphereDataPhysical pot_energy = h*(io_simVars.sim.gravitation*normalization);
		sweet::SphereData_Physical pot_energy = h*h*0.5*normalization;
		sweet::SphereData_Physical kin_energy = h*(u*u+v*v)*(0.5*normalization);

		potential_energy = sphereHelpers_integral.compute_zylinder_integral(pot_energy);
		kinetic_energy = sphereHelpers_integral.compute_zylinder_integral(kin_energy);

		/*
		 * We follow the Williamson et al. equation (137) here
		 */
//		double dummy_energy = compute_zylinder_integral(h*h*(0.5*normalization));
//		io_simVars.diag.total_energy = io_simVars.diag.kinetic_energy + dummy_energy;//io_simVars.diag.potential_energy;

		/*
		 * We follow pot/kin energy here
		 */
		total_energy = kinetic_energy + potential_energy;

		// total vorticity
		// TODO: maybe replace this with the i_vort parameter
		sweet::SphereData_Physical eta(h.sphereDataConfig);
		eta = i_ops->uv_to_vort(u, v).toPhys();

		eta += fg;

		// enstrophy (Williamson paper, equation 138)
		total_potential_enstrophy = 0.5*sphereHelpers_integral.compute_zylinder_integral(eta*eta/h) * normalization;
	}


	void print(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << i_prefix << "DIAGNOSTICS:" << std::endl;
		std::cout << i_prefix << " + total_mass: " << total_mass << std::endl;
		std::cout << i_prefix << " + total_energy: " << total_energy << std::endl;
		std::cout << i_prefix << " + kinetic_energy: " << kinetic_energy << std::endl;
		std::cout << i_prefix << " + potential_energy: " << potential_energy << std::endl;
		std::cout << i_prefix << " + total_potential_enstrophy: " << total_potential_enstrophy << std::endl;
		std::cout << std::endl;
	}

	void printTabularHeader(
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix << "T\tTOTAL_MASS\tPOT_ENERGY\tKIN_ENERGY\tTOT_ENERGY\tPOT_ENSTROPHY\tREL_TOTAL_MASS\tREL_POT_ENERGY\tREL_KIN_ENERGY\tREL_TOT_ENERGY\tREL_POT_ENSTROPHY";
	}

	void printTabularRow(
		double i_current_simulation_time,
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix;

		// Print simulation time, energy and pot enstrophy
		std::cout << i_current_simulation_time << "\t";
		std::cout << total_mass << "\t";
		std::cout << potential_energy << "\t";
		std::cout << kinetic_energy << "\t";
		std::cout << total_energy << "\t";
		std::cout << total_potential_enstrophy << "\t";

		std::cout << (total_mass-ref_total_mass)/total_mass << "\t";
		std::cout << (potential_energy-ref_potential_energy)/potential_energy << "\t";
		std::cout << (kinetic_energy-ref_kinetic_energy)/kinetic_energy << "\t";
		std::cout << (total_energy-total_energy)/total_energy << "\t";
		std::cout << (total_potential_enstrophy-total_potential_enstrophy)/total_potential_enstrophy;
		std::cout << std::endl;
	}
};



#endif /* SRC_PROGRAMS_PDE_SWESPHERE_PDESWESPHEREDIAGNOSTICS_HPP_ */
