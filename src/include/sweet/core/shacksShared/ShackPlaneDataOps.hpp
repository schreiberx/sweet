/*
 *  Created on: Feb 23, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_PLANE_SHACKPLANE_DATA_OPS_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_PLANE_SHACKPLANE_DATA_OPS_HPP_


#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/TransformationPlans.hpp>
#include <sweet/core/StringSplit.hpp>


namespace sweet
{

/**
 * This class stored the discretization-related parameters
 *
 * resolution / timestepping
 */
class ShackPlaneDataOps	:
		public ShackInterface
{
public:
	/**
	 * resolution in physical space (grid cells)
	 */
	int space_res_physical[2] = {0, 0};


	/**
	 * resolution in spectral space (number of modes)
	 */
	int space_res_spectral[2] = {0, 0};


	/**
	 * domain size if running simulation on the plane
	 */
	double plane_domain_size[2] = {1.0, 1.0};


	/**
	 * use spectral differential operators
	 */
	bool space_use_spectral_basis_diffs =
#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_SPHERE_SPECTRAL_SPACE
			true;
#else
			false;
#endif

	/**
	 * Use C-grid staggering
	 */
	bool space_grid_use_c_staggering = false;


	/*
	 * Load / Save plans for SHTNS (useful for reproducibility)
	 *
	 * Can also exist in the SphereDataOps
	 */
	TransformationPlans::TRANSFORMATION_PLAN_CACHE reuse_spectral_transformation_plans = TransformationPlans::QUICK;

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "Discretization:" << std::endl;
		std::cout << "  >Space:" << std::endl;
		std::cout << "	--space-grid-use-c-staggering [0/1]	Use staggering" << std::endl;
		std::cout << "	-N [res]		resolution in x and y direction, default=0" << std::endl;
		std::cout << "	-n [resx]		resolution in x direction, default=0" << std::endl;
		std::cout << "	-m [resy]		resolution in y direction, default=0" << std::endl;
		std::cout << "	-M [modes]		modes in x/y, default=0" << std::endl;
		std::cout << "	-S [0/1]		Control Operator discretization for PlaneData" << std::endl;
		std::cout << "					0: FD, 1: spectral derivatives, default: ";
		std::cout << "	-X [length]	length of simulation domain in x direction, default=1" << std::endl;
		std::cout << "	-Y [width]	width of simulation domain in y direction, default=1" << std::endl;

#if SWEET_USE_PLANE_SPECTRAL_SPACE || SWEET_USE_SPHERE_SPECTRAL_SPACE
std::cout << "1" << std::endl;
#else
std::cout << "0" << std::endl;
#endif

	}

	bool processProgramArguments(ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--space-grid-use-c-staggering", space_grid_use_c_staggering);
		i_pa.getArgumentValueByKey("-S", space_use_spectral_basis_diffs);

		std::string tmp_N;
		if (i_pa.getArgumentValueByKey("-N", tmp_N))
		{

			int c = StringSplit::split2int(tmp_N, &space_res_physical[0], &space_res_physical[1]);
			if (c == 1)
				space_res_physical[1] = space_res_physical[0];
		}

		std::string tmp_M;
		if (i_pa.getArgumentValueByKey("-M", tmp_M))
		{

			int c = StringSplit::split2int(tmp_M, &space_res_spectral[0], &space_res_spectral[1]);
			if (c == 1)
				space_res_spectral[1] = space_res_spectral[0];
		}

		i_pa.getArgumentValueByKey("-X", plane_domain_size[0]);
		i_pa.getArgumentValueByKey("-Y", plane_domain_size[1]);

		/*
		 * We allow multiple --reuse-plans here since this can be also processed by sphere
		 */
		std::string tmp;
		if (i_pa.getArgumentValueByKey("--reuse-plans", tmp, false, false))
		{
			TransformationPlans tp;
			tp.getEnumFromString(tmp, reuse_spectral_transformation_plans);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(tp);
		}

		ERROR_FORWARD_WITH_RETURN_BOOLEAN(i_pa);
	}

	bool validateResolution()
	{
		if (
				(space_res_physical[0] == 0 || space_res_physical[1] == 0)	&&
				(space_res_spectral[0] == 0 || space_res_spectral[1] == 0)
		)
			return error.set("Select physical resolution or spectral modes (use -N (or -n, -m) for physical and -M for spectral");

		return true;
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "DISCRETIZATION:" << std::endl;
		std::cout << i_prefix << " + space_res_physical: " << space_res_physical[0] << " x " << space_res_physical[1] << std::endl;
		std::cout << i_prefix << " + space_res_spectral: " << space_res_spectral[0] << " x " << space_res_spectral[1] << std::endl;
		std::cout << i_prefix << " + space_use_spectral_basis_diffs: " << space_use_spectral_basis_diffs << std::endl;
		std::cout << i_prefix << " + space_grid_use_c_staggering: " << space_grid_use_c_staggering << std::endl;
		std::cout << i_prefix << " + plane_dealiasing (compile time): " <<
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		1
#else
		0
#endif
		<< std::endl;
		std::cout << i_prefix << " + domain_size: " << plane_domain_size[0] << " x " << plane_domain_size[1] << std::endl;
		std::string tmp;
		TransformationPlans tp;
		tp.getStringFromEnum(reuse_spectral_transformation_plans, tmp);
		std::cout << i_prefix << " + reuse_spectral_transformation_plans: " << tmp << std::endl;
		std::cout << i_prefix << std::endl;
	}
};

}

#endif
