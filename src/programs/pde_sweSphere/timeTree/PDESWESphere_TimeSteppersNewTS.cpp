#include "PDESWESphere_TimeSteppersNewTS.hpp"


#include "PDESWESphere_ln.hpp"

#include "PDESWESphere_l.hpp"
#include "PDESWESphere_lg.hpp"
#include "PDESWESphere_lc.hpp"

#include "PDESWESphere_n.hpp"
#include "PDESWESphere_na_uv.hpp"
#include "PDESWESphere_nr_uv.hpp"
#include "PDESWESphere_na_vd.hpp"
#include "PDESWESphere_nr_vd.hpp"


bool PDESWESphere_TimeSteppersNewTS::setup_1_registerAllTimesteppers()
{

	/*
	 * Registration of all possible PDE terms
	 */
	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_ln>();

	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_l>();
	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_lg>();
	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_lc>();

	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_n>();
	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_na_uv>();
	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_nr_uv>();
	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_na_vd>();
	pdeTerm_registry.registerTimeTreeNode<PDESWESphere_nr_vd>();

	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeTerm_registry);

	/*
	 * Register time steppers
	 */
	sweet::DESolver_TimeStepperRegistryAll registryAll;
	registryAll.registerAll(timeStepper_registry);

	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(registryAll);

	return true;
}


PDESWESphere_TimeSteppersNewTS::PDESWESphere_TimeSteppersNewTS()
{
}

bool PDESWESphere_TimeSteppersNewTS::setup_2_timestepper(
		const std::string &i_timestepping_method,
		sweet::ShackDictionary *i_shackDict,
		sweet::SphereOperators *io_ops,
		const PDESWESphere_DataContainer &i_U
)
{
	/*
	 * Setup time stepping string parser and parse it
	 */
	sweet::DESolver_TimeSteppingStringParser tsStringParser;
	tsStringParser.genTimeSteppingTree(
			i_timestepping_method,
			tsTree
		);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tsStringParser);
	tsTree.print();

	/*
	 * Ready to assemble time stepper
	 */
	sweet::DESolver_TimeStepping_Assemblation tssa;
	tssa.setup(pdeTerm_registry, timeStepper_registry);
	tssa.assembleTimeStepperByTree(
		tsTree,
		timeIntegrator
	);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tssa);

	timeIntegrator->shackRegistration(i_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeIntegrator);

	deSolver_Config.myDataContainer = &i_U;
	deSolver_Config.ops = io_ops;
	timeIntegrator->setupConfig(deSolver_Config);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeIntegrator);

	/*
	 * Set time step size
	 */
	timeIntegrator->setTimeStepSize(0.1);

	return true;
}

void PDESWESphere_TimeSteppersNewTS::clear()
{
}


PDESWESphere_TimeSteppersNewTS::~PDESWESphere_TimeSteppersNewTS()
{
}
