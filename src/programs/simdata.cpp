/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * disabled_MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/
 * disabled_MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/time
 * disabled_MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/benchmarks
 *
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 * disabled_MULE_SCONS_OPTIONS: --fortran-source=enable
 * disabled_MULE_SCONS_OPTIONS: --lapack=enable
 */

#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/ErrorBase.hpp>

#include "simdata/MyDataContainer.hpp"

#include "simdata/PDESolver_PDETerm_Registry.hpp"
#include "simdata/MyPDETerm_lg.hpp"
#include "simdata/MyPDETerm_lc.hpp"

#include "simdata/PDESolver_TimeStepper_Base.hpp"
#include "simdata/PDESolver_TimeStepper_Registry.hpp"
#include "simdata/MyTimeStepper_ExplicitRungeKutta.hpp"

#include "simdata/PDESolver_TimeStepping_StringParser.hpp"
#include "simdata/PDESolver_TimeStepping_Tree.hpp"
#include "simdata/PDESolver_TimeStepping_Assemblation.hpp"


class ShackTimeSteppingMethod :
		public sweet::ShackInterface
{
public:
	/// String of time stepping method
	std::string timestepping_method;

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "Time Discretization:" << std::endl;
		std::cout << "	--timestepping-method [string]	String of time stepping method" << std::endl;

	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--timestepping-method", timestepping_method);

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
		std::cout << std::endl;
	}
};



class PDESolver
{
public:
	sweet::ErrorBase error;

private:
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	ShackTimeSteppingMethod *shackTimeSteppingMethod;
	sweet::SphereData_Config sphereDataConfig;
	sweet::SphereOperators ops;


	MyDataContainer U;
	MyDataContainer U_tmp;

public:
	PDESolver()	:
		shackSphereDataOps(nullptr)
	{
	}

	bool setup_1_shacks(int argc, char *argv[])
	{
		// setup shacks
		shackProgArgDict.setup(argc, argv);
		shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
		shackTimeSteppingMethod = shackProgArgDict.getAutoRegistration<ShackTimeSteppingMethod>();
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	bool setup_2_config()
	{
		// setup sphere data config
		sphereDataConfig.setupAuto(shackSphereDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereDataConfig);

		// setup sphere operators
		ops.setup(&sphereDataConfig, shackSphereDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

		return true;
	}

	bool setup_3_data()
	{
		U.setup(&sphereDataConfig);
		U.phi.spectral_set_zero();
		U.vrt.spectral_set_zero();
		U.div.spectral_set_zero();

		U_tmp.setup(&sphereDataConfig);

		return true;
	}

	std::shared_ptr<sweet::PDESolver_TimeStepper_Base> timeStepper;

	bool setup_4_timestepper(int rk_order)
	{
#if 1
		/*
		 * Setup time stepping string parser and parse it
		 */
		sweet::PDESolver_TimeSteppingStringParser tsStringParser;
		sweet::PDESolver_TimeStepping_Tree tsTree;

		tsStringParser.genTimeSteppingTree(
				shackTimeSteppingMethod->timestepping_method,
				tsTree
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tsStringParser);
		tsTree.print();

		/*
		 * Register PDE Terms
		 */
		sweet::PDESolver_PDETerm_Registry pdeTerm_registry;
		pdeTerm_registry.registerPDETerm<MyPDETerm_lg>();
		pdeTerm_registry.registerPDETerm<MyPDETerm_lc>();

		/*
		 * Register time steppers
		 */
		sweet::PDESolver_TimeStepper_Registry timeStepper_registry;
		timeStepper_registry.registerTimeStepper<MyTimeStepper_ExplicitRungeKutta>();

		/*
		 * Ready to assemble time stepper
		 */
		sweet::PDESolver_TimeStepping_Assemblation tssa;
		tssa.setup(pdeTerm_registry, timeStepper_registry);
		tssa.assembleTimeStepperByTree(
			tsTree,
			timeStepper
		);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tssa);

#else

		std::shared_ptr<sweet::PDESolver_PDETerm_Base> pde_term_lg;

		/*
		 * Register
		 */
		pde_term_lg = std::shared_ptr<sweet::PDESolver_PDETerm_Base>(new MyPDETerm_lg);
		pde_term_lg->shackRegistration(&shackProgArgDict);
		pde_term_lg->setup(&ops, U);

		timeStepper = std::make_shared<MyTimeStepper_ExplicitRungeKutta>();
		timeStepper->shackRegistration(&shackProgArgDict);

		MyTimeStepper_ExplicitRungeKutta *_ERK = static_cast<MyTimeStepper_ExplicitRungeKutta*>(timeStepper.get());
		_ERK->setupWithArguments(pde_term_lg, rk_order);
		timeStepper->setTimestepSize(0.1);
		timeStepper->setupDataContainers(U);

#endif

		return true;
	}

	void clear()
	{
	}

	void doTimestep()
	{
		timeStepper->eval_timeIntegration(U, U_tmp, 0);
		U.swap(U_tmp);
	}
};


bool runTests()
{
	sweet::PDESolver_TimeSteppingStringParser tssParser;
	sweet::PDESolver_TimeStepping_Tree tsTree;

	{
		std::string tmp = "foo()";
		std::cout << "***" << std::endl;
		std::cout << tmp << std::endl;
		std::cout << "***" << std::endl;
		tssParser.genTimeSteppingTree(tmp, tsTree);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(tssParser);
		tsTree.print();
		std::cout << "***" << std::endl;
	}

	{
		std::string tmp = "foo(1, 2, 3)";
		std::cout << "***" << std::endl;
		std::cout << tmp << std::endl;
		std::cout << "***" << std::endl;
		tssParser.genTimeSteppingTree(tmp, tsTree);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(tssParser);
		tsTree.print();
		std::cout << "***" << std::endl;
	}

	{
		std::string tmp = "foo(a=1, b=3, c=foo)";
		std::cout << "***" << std::endl;
		std::cout << tmp << std::endl;
		std::cout << "***" << std::endl;
		tssParser.genTimeSteppingTree(tmp, tsTree);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(tssParser);
		tsTree.print();
		std::cout << "***" << std::endl;
	}

	{
		std::string tmp = "foo(bar(), fuzz(a=b))";
		std::cout << "***" << std::endl;
		std::cout << tmp << std::endl;
		std::cout << "***" << std::endl;
		tssParser.genTimeSteppingTree(tmp, tsTree);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(tssParser);
		tsTree.print();
		std::cout << "***" << std::endl;
	}

	{
		std::string tmp = "foo(a=bar(a=b), barx())";
		std::cout << "***" << std::endl;
		std::cout << tmp << std::endl;
		std::cout << "***" << std::endl;
		tssParser.genTimeSteppingTree(tmp, tsTree);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(tssParser);
		tsTree.print();
		std::cout << "***" << std::endl;
	}

	{
		std::string tmp = "FunctionNameAnd123SomeNumbers(these, are=parameters, and, can, also, be, a, recursive, function(yadda))";
		std::cout << "***" << std::endl;
		std::cout << tmp << std::endl;
		std::cout << "***" << std::endl;
		tssParser.genTimeSteppingTree(tmp, tsTree);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(tssParser);
		tsTree.print();
		std::cout << "***" << std::endl;
	}

	return true;
}

int main(int argc, char *argv[])
{
	if (!runTests())
		return -1;

	for (int rk_order = 1; rk_order <= 4; rk_order++)
	{
		std::cout << "testing rk_order: " << rk_order << std::endl;
		PDESolver pdeSolver;

		pdeSolver.setup_1_shacks(argc, argv);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(pdeSolver);

		pdeSolver.setup_2_config();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(pdeSolver);

		pdeSolver.setup_3_data();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(pdeSolver);

		pdeSolver.setup_4_timestepper(rk_order);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(pdeSolver);

		pdeSolver.doTimestep();
	}

	return 0;
}
