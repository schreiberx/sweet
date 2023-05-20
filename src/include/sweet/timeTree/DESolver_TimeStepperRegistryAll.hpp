/*
 * DESolver_TimeStepperRegistryAll.hpp
 *
 *  Created on: May 18, 2023
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_TIMESTEPPERREGISTRY_HPP_
#define SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_TIMESTEPPERREGISTRY_HPP_


#include <sweet/core/ErrorBase.hpp>
#include "DESolver_TimeTreeNode_Registry.hpp"

#include <sweet/timeTree/DESolver_TimeStepper_AddDETerms.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_NegDETerms.hpp>

#include <sweet/timeTree/DESolver_TimeStepper_ExplicitRungeKutta.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_ImplicitRungeKutta.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_StrangSplitting.hpp>


namespace sweet
{

/*
 * This is a class to register all time steppers
 */
class DESolver_TimeStepperRegistryAll
{
public:
	ErrorBase error;

public:
	DESolver_TimeStepperRegistryAll()
	{
	}

	bool registerAll(
			sweet::DESolver_TimeTreeNode_Registry &o_timeStepper_registry
	)
	{
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_AddDETerms>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_NegDETerms>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_ExplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_ImplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_StrangSplitting>();

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(o_timeStepper_registry);
		return true;
	}
};

}
#endif /* SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_TIMESTEPPERREGISTRY_HPP_ */
