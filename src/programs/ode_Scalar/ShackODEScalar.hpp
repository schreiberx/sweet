/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_COMMON_SHACK_ODESCALARPARAMETERSCOMMON_HPP_
#define SRC_PROGRAMS_SWE_COMMON_SHACK_ODESCALARPARAMETERSCOMMON_HPP_

#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/ErrorBase.hpp>


class ShackODEScalar	:
		public sweet::ShackInterface
{
public:

	//////*
	///// * Initial solution
	///// */
	/////double u0;

	//////*
	///// * Parameters a and b from u_t = a*sin(u) = b*sin(t)
	///// */
	/////double param_a;
	/////double param_b;


	/*
	 * Compute errors compared to analytical solution
	 */
	bool compute_errors = false;


	/*
	 * Compute diagnostics
	 */
	bool compute_diagnostics = false;


	/*
	 * Check for instabilities and stop
	 */
	bool instability_checks = false;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "ODEScalarParameters parameters:" << std::endl;
		////std::cout << i_prefix << "	--u0=[float]			Initial solution" << std::endl;
		////std::cout << i_prefix << "	--param-a=[float]		ODE parameter a" << std::endl;
		////std::cout << i_prefix << "	--param-b=[float]		ODE parameter b" << std::endl;
		std::cout << i_prefix << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << i_prefix << "	--compute-diagnostics [bool]	Compute diagnostics" << std::endl;
		std::cout << i_prefix << "	--instability-checks=[bool]			Check for instabilities (default: 0)" << std::endl;
	}


	void printShack(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "ODEScalar parameters:" << std::endl;
		////std::cout << i_prefix << " + u0: " << u0 << std::endl;
		////std::cout << i_prefix << " + param_a: " << param_a << std::endl;
		////std::cout << i_prefix << " + param_b: " << param_b << std::endl;
		std::cout << i_prefix << " + compute_errors: " << compute_errors << std::endl;
		std::cout << i_prefix << " + compute_diagnostics: " << compute_diagnostics << std::endl;
		std::cout << i_prefix << " + instability_checks: " << instability_checks << std::endl;
	}



	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		////i_pa.getArgumentValueByKey("-u0", u0);
		////i_pa.getArgumentValueByKey("-param-a", param_a);
		////i_pa.getArgumentValueByKey("-param-b", param_b);

		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);
		i_pa.getArgumentValueByKey("--compute-diagnostics", compute_diagnostics);
		i_pa.getArgumentValueByKey("--instability-checks", instability_checks);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_pa);
		return true;
	}
};



#endif
