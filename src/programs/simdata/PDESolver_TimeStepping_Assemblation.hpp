#ifndef SRC_INCLUDE_SWEET_PDESOLVER_TIMESTEPPINGSTRING_ASSEMBLATION_HPP_
#define SRC_INCLUDE_SWEET_PDESOLVER_TIMESTEPPINGSTRING_ASSEMBLATION_HPP_

#include <string>
#include <vector>
#include <ostream>
#include <memory>

#include <sweet/core/ErrorBase.hpp>
#include "PDESolver_TimeStepping_Tree.hpp"
#include "PDESolver_TimeStepper_Registry.hpp"


// Forward declaration of Base
// It's included at the end of this file
namespace sweet {
	class PDESolver_TimeStepper_Base;
}

namespace sweet
{

/*
 * Assembles a time stepper from a time stepping string
 */
class PDESolver_TimeStepping_Assemblation
{
public:
	ErrorBase error;

private:
	PDESolver_PDETerm_Registry *pdeTermsRegistry;

	PDESolver_TimeStepper_Registry *timeSteppersRegistry;

public:
	PDESolver_TimeStepping_Assemblation()	:
		pdeTermsRegistry(nullptr),
		timeSteppersRegistry(nullptr)
	{
	}


	void clear()
	{
	}


	bool setup(
		PDESolver_PDETerm_Registry &i_pdeTerms,
		PDESolver_TimeStepper_Registry &i_timeSteppers
	)
	{
		pdeTermsRegistry = &i_pdeTerms;
		timeSteppersRegistry = &i_timeSteppers;

		return true;
	}


public:
	/*
	 * Setup the full time stepping tree
	 */
	bool assembleTimeStepperByTree(
			PDESolver_TimeStepping_Tree &i_tree,
			std::shared_ptr<PDESolver_TimeStepper_Base> &o_timestepper
	)
	{
		if (pdeTermsRegistry == nullptr || timeSteppersRegistry == nullptr)
			return error.set("You need to call setup(...) before the assemblation");

		return assembleTimeStepperByFunction(
				i_tree.mainFunction,
				o_timestepper
			);
	}


public:
	/*
	 * Setup the time stepper for a given function and return it
	 */
	bool assembleTimeStepperByFunction(
		std::shared_ptr<PDESolver_TimeStepping_Tree::Function> &i_function,
		std::shared_ptr<PDESolver_TimeStepper_Base> &o_timestepper
	)
	{
		/*
		 * Step 1) Search for implementation of time stepper of this particular function
		 */
		timeSteppersRegistry->getTimeStepperNewInstance(
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

	bool assemblePDETermByString(
			const std::string i_pde_string,
			std::shared_ptr<PDESolver_PDETerm_Base> &o_pdeTerm
	)
	{
		pdeTermsRegistry->getPDETermInstance(i_pde_string, o_pdeTerm);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*pdeTermsRegistry);

		return true;
	}

#if 0
	friend
	std::ostream&
	operator<<(std::ostream &io_os, const PDESolver_TimeSteppingStringParser &i_pa)
	{
		return io_os;
	}
#endif
};

}

// Included here due to forward declaration
#include "PDESolver_TimeStepper_Base.hpp"

#endif
