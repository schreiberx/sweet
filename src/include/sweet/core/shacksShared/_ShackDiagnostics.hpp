/*
 * ShackDiagnostics.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKDIAGNOSTICS_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKDIAGNOSTICS_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>



/**
 * Diagnostic variables
 */
class ShackDiagnostics	:
	public ShackInterface
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

	bool processProgramArguments(ProgramArguments &i_pa)
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
};




#endif
