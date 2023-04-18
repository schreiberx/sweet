/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_SWEET_SHACK_PDE_SWE_PLANE_MORI_ZWANZIG_TIME_DISCRETIZATION_HPP_
#define SRC_INCLUDE_SWEET_SHACK_PDE_SWE_PLANE_MORI_ZWANZIG_TIME_DISCRETIZATION_HPP_


#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/StringSplit.hpp>


class ShackPDESWEPlaneMoriZwanzigTimeDiscretization	:
		public sweet::ShackInterface
{
public:

	/// String of time stepping method
	/// See doc/swe/swe_plane_timesteppings
	std::string timestepping_method;
	std::string timestepping_method_P;
	std::string timestepping_method_Q;

	/// Order of time stepping
	std::string timestepping_order;
	int timestepping_order_P = -1;
	int timestepping_order_Q = -1;

	/// Order of 2nd time stepping which might be used
	std::string timestepping_order2;
	int timestepping_order2_P = -1;
	int timestepping_order2_Q = -1;

	///// Number of iterations for semi-Lagrangian methods
	//std::string semi_lagrangian_max_iterations;
	//int semi_lagrangian_max_iterations_SP = 2;
	//int semi_lagrangian_max_iterations_SQ = 2;
	//int semi_lagrangian_max_iterations_FQ = 2;

	///// Convergence threshold for semi-Lagrangian methods (set to -1 to ignore error)
	//std::string semi_lagrangian_convergence_threshold;
	//double semi_lagrangian_convergence_threshold_SP = -1;
	//double semi_lagrangian_convergence_threshold_SQ = -1;
	//double semi_lagrangian_convergence_threshold_FQ = -1;

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "Time Discretization:" << std::endl;
		std::cout << "	--timestepping-method [string]	String of time stepping method" << std::endl;
		std::cout << "	--timestepping-order [int]			Specify the order of the time stepping" << std::endl;
		std::cout << "	--timestepping-order2 [int]			Specify the order of the time stepping" << std::endl;
		//std::cout << "	--semi-lagrangian-max-iterations [int]		Number of max. iterations during semi-Lagrangian time integration" << std::endl;
		//std::cout << "	--semi-lagrangian-convergence-threshold [float]	Threshold to stop iterating, Use -1 to disable" << std::endl;

	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		///i_pa.getArgumentValueByKey("--timestepping-method", timestepping_method);
		///i_pa.getArgumentValueByKey("-R", timestepping_order);
		///i_pa.getArgumentValueByKey("--timestepping-order", timestepping_order);
		//i_pa.getArgumentValueByKey("--timestepping-order2", timestepping_order2);
		//i_pa.getArgumentValueByKey("--semi-lagrangian-max-iterations", semi_lagrangian_max_iterations);
		//i_pa.getArgumentValueByKey("--semi-lagrangian-convergence-threshold", semi_lagrangian_convergence_threshold);

		i_pa.getArgumentValueByKey("--MZ-timestepping-method-P", timestepping_method_P);
		i_pa.getArgumentValueByKey("--MZ-timestepping-method-Q", timestepping_method_Q);
		i_pa.getArgumentValueByKey("--MZ-timestepping-order-P", timestepping_order_P);
		i_pa.getArgumentValueByKey("--MZ-timestepping-order-Q", timestepping_order_Q);
		i_pa.getArgumentValueByKey("--MZ-timestepping-order2-P", timestepping_order2_P);
		i_pa.getArgumentValueByKey("--MZ-timestepping-order2-Q", timestepping_order2_Q);


		////std::vector<std::string> split;

		////split = StringSplit::split(timestepping_method, ",");
		////if (split.size() == 1)
		////{
		////	std::string v = split[0];
		////	timestepping_method_P = v;
		////	timestepping_method_Q = v;
		////}
		////else if (split.size() == 2)
		////{
		////	timestepping_method_P = split[0];
		////	timestepping_method_Q = split[1];
		////}
		////else
		////	SWEETError("Invalid number of arguments in timestepping_method");

		/////split = StringSplit::split(timestepping_order, ",");
		/////if (split.size() == 1)
		/////{
		/////	int v = std::stoi(split[0]);
		/////	timestepping_order_P = v;
		/////	timestepping_order_Q = v;
		/////}
		/////else if (split.size() == 2)
		/////{
		/////	timestepping_order_P = std::stoi(split[0]);
		/////	timestepping_order_Q = std::stoi(split[1]);
		/////}
		/////else
		/////	SWEETError("Invalid number of arguments in timestepping_order");

		/////split = StringSplit::split(timestepping_order2, ",");
		/////if (split.size() == 1)
		/////{
		/////	int v = std::stoi(split[0]);
		/////	timestepping_order2_P = v;
		/////	timestepping_order2_Q = v;
		/////}
		/////else if (split.size() == 2)
		/////{
		/////	timestepping_order2_P = std::stoi(split[0]);
		/////	timestepping_order2_Q = std::stoi(split[1]);
		/////}
		/////else
		/////	SWEETError("Invalid number of arguments in timestepping_order2");

		////split = StringSplit::split(semi_lagrangian_max_iterations, ",");
		////if (split.size() == 1)
		////{
		////	int v = std::stoi(split[0]);
		////	semi_lagrangian_max_iterations_SP = v;
		////	semi_lagrangian_max_iterations_SQ = v;
		////	semi_lagrangian_max_iterations_FQ = v;
		////}
		////else if (split.size() == 3)
		////{
		////	semi_lagrangian_max_iterations_SP = std::stoi(split[0]);
		////	semi_lagrangian_max_iterations_SQ = std::stoi(split[1]);
		////	semi_lagrangian_max_iterations_FQ = std::stoi(split[2]);
		////}
		////else
		////	SweetError("Invalid number of arguments in semo_lagrangian_max_iterations");

		////split = StringSplit::split(semi_lagrangian_convergence_threshold, ",");
		////if (split.size() == 1)
		////{
		////	int v = std::stod(split[0]);
		////	semi_lagrangian_convergence_threshold_SP = v;
		////	semi_lagrangian_convergence_threshold_SQ = v;
		////	semi_lagrangian_convergence_threshold_FQ = v;
		////}
		////else if (split.size() == 3)
		////{
		////	semi_lagrangian_convergence_threshold_SP = std::stod(split[0]);
		////	semi_lagrangian_convergence_threshold_SQ = std::stod(split[1]);
		////	semi_lagrangian_convergence_threshold_FQ = std::stod(split[2]);
		////}
		////else
		////	SweetError("Invalid number of arguments in semo_lagrangian_convergence_threshold");


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
		////std::cout << " + semi_lagrangian_max_iterations: " << semi_lagrangian_max_iterations << std::endl;
		////std::cout << " + semi_lagrangian_convergence_threshold: " << semi_lagrangian_convergence_threshold << std::endl;
		std::cout << std::endl;
	}
};




#endif
