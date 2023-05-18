

/*
 * Always include these classes due to forward delaration stuff
 */

#ifndef SRC_INCLUDE_SWEET_DESOLVER_TIMESTEPPINGSTRING_ASSEMBLATION_HPP_
#define SRC_INCLUDE_SWEET_DESOLVER_TIMESTEPPINGSTRING_ASSEMBLATION_HPP_

#include <string>
#include <vector>
#include <ostream>
#include <memory>
#include <sweet/core/ErrorBase.hpp>
#include "DESolver_TimeStepping_Tree.hpp"
#include "DESolver_TimeTreeNode_Registry.hpp"



// Forward declaration of Base
// It's included at the end of this file
namespace sweet {
	class DESolver_TimeTreeNode_Base;
}

namespace sweet
{

/*
 * Assembles a time stepper from a time stepping string
 */
class DESolver_TimeStepping_Assemblation
{
public:
	ErrorBase error;

private:
	DESolver_TimeTreeNode_Registry *deTermsRegistry;

	DESolver_TimeTreeNode_Registry *timeSteppersRegistry;

public:
	DESolver_TimeStepping_Assemblation()	:
		deTermsRegistry(nullptr),
		timeSteppersRegistry(nullptr)
	{
	}


	void clear()
	{
	}


	bool setup(
		DESolver_TimeTreeNode_Registry &i_pdeTerms,
		DESolver_TimeTreeNode_Registry &i_timeSteppers
	)
	{
		deTermsRegistry = &i_pdeTerms;
		timeSteppersRegistry = &i_timeSteppers;

		return true;
	}


public:
	/*
	 * Setup the full time stepping tree
	 */
	bool assembleTimeStepperByTree(
			DESolver_TimeStepping_Tree &i_tree,
			std::shared_ptr<DESolver_TimeTreeNode_Base> &o_timestepper
	)
	{
		if (deTermsRegistry == nullptr || timeSteppersRegistry == nullptr)
			return error.set("You need to call setup(...) before the assemblation");

		return assembleTimeTreeNodeByFunction(
				i_tree.mainFunction,
				o_timestepper
			);
	}


public:
	/*
	 * Setup the time stepper for a given function and return it
	 */
	bool assembleTimeTreeNodeByFunction(
		std::shared_ptr<DESolver_TimeStepping_Tree::Function> &i_function,
		std::shared_ptr<DESolver_TimeTreeNode_Base> &o_timestepper
	)
	{
		/*
		 * Step 1) Search for implementation of time stepper of this particular function
		 */
		timeSteppersRegistry->getTimeTreeNodeNewInstance(
				i_function->function_name,
				o_timestepper
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeSteppersRegistry);

		/*
		 * Step 2) Call the setup routine of the time stepper because
		 * Only the time steppers knows about how to process its arguments.
		 *
		 * We also hand over this class since there could be other
		 *
		 */
		o_timestepper->setupFunction(i_function, *this);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*(o_timestepper.get()));
		return true;
	}


public:
	/*
	 * Setup the time stepper for a given function and return it
	 */
	bool assembleTimeTreeNodeByName(
		const std::string i_nodeName,
		std::shared_ptr<DESolver_TimeTreeNode_Base> &o_timestepper
	)
	{
		/*
		 * Step 1) Search for implementation of time stepper of this particular function
		 */

		// we first search for the de terms
		if (!deTermsRegistry->getTimeTreeNodeNewInstance(
				i_nodeName,
				o_timestepper
			)
		)
		{
			deTermsRegistry->error.reset();

			timeSteppersRegistry->getTimeTreeNodeNewInstance(
					i_nodeName,
					o_timestepper
				);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timeSteppersRegistry);
		}

		return true;
	}

#if 0
	bool assembleDETermByString(
			const std::string i_pde_string,
			std::shared_ptr<DESolver_TimeTreeNode_Base> &o_pdeTerm
	)
	{
		deTermsRegistry->getTimeTreeNodeNewInstance(i_pde_string, o_pdeTerm);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*deTermsRegistry);

		return true;
	}
#endif

#if 0
	friend
	std::ostream&
	operator<<(std::ostream &io_os, const DESolver_TimeSteppingStringParser &i_pa)
	{
		return io_os;
	}
#endif
};

}

// Included here due to forward declaration
#include "DESolver_TimeTreeNode_Base.hpp"


#endif
