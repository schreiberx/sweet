/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_PLANE_SHACK_SPHERE_DATA_OPS_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_PLANE_SHACK_SPHERE_DATA_OPS_HPP_


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
class ShackSphereDataOps	:
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
	 * Earth radius for simulations on the sphere
	 */
	double sphere_radius = 6.37122e6;


	/*
	 * Verbosity level for SH setup
	 */
	int sh_setup_verbosity = 0;

	/*
	 * Number of threads during setup (useful for nested parallelization)
	 */
	int sh_setup_num_threads = -1;



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
		std::cout << "	-N [res]		resolution in x and y direction, default=0" << std::endl;
		std::cout << "	-n [resx]		resolution in x direction, default=0" << std::endl;
		std::cout << "	-m [resy]		resolution in y direction, default=0" << std::endl;
		std::cout << "	-M [modes]		modes in x/y, default=0" << std::endl;
		std::cout << "	-X [length]	length of simulation domain in x direction, default=1" << std::endl;
		std::cout << "	-Y [width]	width of simulation domain in y direction, default=1" << std::endl;
		std::cout << "	--sh-setup-num-threads [int]	Number of threads for SH setup, default=-1 (disabled)" << std::endl;
		std::cout << "	--sh-setup-verbosity [int]	Verbosity during SH setup, default=0" << std::endl;
	}

	bool processProgramArguments(ProgramArguments &i_pa)
	{
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

		/*
		 * Sphere
		 */
		i_pa.getArgumentValueByKey("-a", sphere_radius);


		i_pa.getArgumentValueByKey("--sh-setup-num-threads", sh_setup_num_threads);
		i_pa.getArgumentValueByKey("--sh-setup-verbosity", sh_setup_verbosity);

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

		std::string tmp;
		TransformationPlans tp;
		tp.getStringFromEnum(reuse_spectral_transformation_plans, tmp);
		std::cout << i_prefix << " + reuse_spectral_transformation_plans: " << tmp << std::endl;
		std::cout << i_prefix << " + sphere_radius: " << sphere_radius << std::endl;
		std::cout << i_prefix << " + sh_setup_num_threads: " << sh_setup_num_threads << std::endl;
		std::cout << i_prefix << " + sh_setup_verbosity: " << sh_setup_verbosity << std::endl;

		std::cout << i_prefix << std::endl;
	}
};

}

#endif
