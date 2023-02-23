/*
 *  Created on: Feb 23, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACK_TIME_DISCRETIZATION_HPP_
#define SRC_INCLUDE_SWEET_SHACK_TIME_DISCRETIZATION_HPP_


#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>


/**
 * This class stored the discretization-related parameters
 *
 * resolution / timestepping
 */
class ShackTimeDiscretization	:
		public sweet::ShackInterface
{
public:
	/// Leapfrog: Robert Asselin filter
	double timestepping_leapfrog_robert_asselin_filter = 0;

	/// Crank-Nicolson filter
	double timestepping_crank_nicolson_filter = 0.5;

	/// Number of iterations for semi-Lagrangian methods
	int semi_lagrangian_max_iterations = 2;

	/// Method to use for computing depature points
	std::string semi_lagrangian_departure_point_method = "settls";

	/// Use limiter for higher-order SL interpolation (cubic)
	bool semi_lagrangian_interpolation_limiter = false;

	/// Create pseudo points at poles for interpolation
	bool semi_lagrangian_sampler_use_pole_pseudo_points = false;

	/// Convergence threshold for semi-Lagrangian methods (set to -1 to ignore error)
	double semi_lagrangian_convergence_threshold = -1;

	/// Use accurate spherical geometry (???) or approximation (Ritchie 1995)
	double semi_lagrangian_approximate_sphere_geometry = 0;


	/// String of time stepping method
	/// See doc/swe/swe_plane_timesteppings
	std::string timestepping_method;

	/// Order of time stepping
	int timestepping_order = -1;

	/// Order of 2nd time stepping which might be used
	int timestepping_order2 = -1;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "Time Discretization:" << std::endl;
		std::cout << "	--timestepping-method [string]	String of time stepping method" << std::endl;
		std::cout << "	--timestepping-order [int]			Specify the order of the time stepping" << std::endl;
		std::cout << "	--timestepping-order2 [int]			Specify the order of the time stepping" << std::endl;
		std::cout << "	--leapfrog-robert-asselin-filter [0;1]		Damping parameter for Robert-Asselin filter" << std::endl;
		std::cout << "	--normal-mode-analysis-generation [0;1;2;3]	Generate output data for normal mode analysis" << std::endl;
		std::cout << "							0: don't generate" << std::endl;
		std::cout << "							1: generate in physical space" << std::endl;
		std::cout << "							2: generate in spectral space" << std::endl;
		std::cout << "							3: generate in spectral space with complex matrix" << std::endl;
		std::cout << "	--semi-lagrangian-max-iterations [int]		Number of max. iterations during semi-Lagrangian time integration" << std::endl;
		std::cout << "	--semi-lagrangian-departure-point-method [str]		'settls' (default), 'midpoint', 'std'" << std::endl;
		std::cout << "	--semi-lagrangian-sampler-use-pole-pseudo-points [bool]" << std::endl;
		std::cout << "								Use pseudo poles for sampling" << std::endl;
		std::cout << "									false (default)" << std::endl;
		std::cout << "	--semi-lagrangian-interpolation-limiter [bool]	Use limiter for cubic interpolation" << std::endl;
		std::cout << "	--semi-lagrangian-convergence-threshold [float]	Threshold to stop iterating, Use -1 to disable" << std::endl;
		std::cout << "	--semi-lagrangian-approximate-sphere-geometry [int]	0: no approximation, 1: Richies approximation, default: 0" << std::endl;

	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--timestepping-method", timestepping_method);
		i_pa.getArgumentValueByKey("-R", timestepping_order);
		i_pa.getArgumentValueByKey("--timestepping-order", timestepping_order);
		i_pa.getArgumentValueByKey("--timestepping-order2", timestepping_order2);
		i_pa.getArgumentValueByKey("--leapfrog-robert-asselin-filter", timestepping_leapfrog_robert_asselin_filter);
		i_pa.getArgumentValueByKey("--crank-nicolson-filter", timestepping_crank_nicolson_filter);
		i_pa.getArgumentValueByKey("--semi-lagrangian-max-iterations", semi_lagrangian_max_iterations);
		i_pa.getArgumentValueByKey("--semi-lagrangian-departure-point-method", semi_lagrangian_departure_point_method);
		i_pa.getArgumentValueByKey("--semi-lagrangian-sampler-use-pole-pseudo-points", semi_lagrangian_sampler_use_pole_pseudo_points);
		i_pa.getArgumentValueByKey("--semi-lagrangian-interpolation-limiter", semi_lagrangian_interpolation_limiter);
		i_pa.getArgumentValueByKey("--semi-lagrangian-convergence-threshold", semi_lagrangian_convergence_threshold);
		i_pa.getArgumentValueByKey("--semi-lagrangian-approximate-sphere-geometry", semi_lagrangian_approximate_sphere_geometry);


		if (i_pa.error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		return true;

	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "TIME DISCRETIZATION:" << std::endl;
		std::cout << " + timestepping_method: " << timestepping_method << std::endl;
		std::cout << " + timestepping_order: " << timestepping_order << std::endl;
		std::cout << " + timestepping_order2: " << timestepping_order2 << std::endl;
		std::cout << " + timestepping_leapfrog_robert_asselin_filter: " << timestepping_leapfrog_robert_asselin_filter << std::endl;
		std::cout << " + timestepping_crank_nicolson_filter: " << timestepping_crank_nicolson_filter << std::endl;
		std::cout << " + semi_lagrangian_max_iterations: " << semi_lagrangian_max_iterations << std::endl;
		std::cout << " + semi_lagrangian_departure_point_method: " << semi_lagrangian_departure_point_method << std::endl;
		std::cout << " + semi_lagrangian_interpolation_limiter: " << semi_lagrangian_interpolation_limiter << std::endl;
		std::cout << " + semi_lagrangian_sampler_use_pole_pseudo_points: " << semi_lagrangian_sampler_use_pole_pseudo_points << std::endl;
		std::cout << " + semi_lagrangian_convergence_threshold: " << semi_lagrangian_convergence_threshold << std::endl;
		std::cout << " + semi_lagrangian_approximate_sphere_geometry: " << semi_lagrangian_approximate_sphere_geometry << std::endl;
		std::cout << " + plane_dealiasing (compile time): " <<
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		1
#else
		0
#endif
		<< std::endl;
		std::cout << std::endl;
	}
};




#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKDISCRETIZATION_HPP_ */
