/*
 * DESolver_TimeStepperRegistryAll.hpp
 *
 *  Created on: May 18, 2023
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_TIMESTEPPERREGISTRY_HPP_
#define SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_TIMESTEPPERREGISTRY_HPP_


#include "DESolver_TimeTreeNode_Registry.hpp"

#include <sweet/timeTree/DESolver_TimeStepper_ExplicitRungeKutta.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_ImplicitRungeKutta.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_AddDETerms.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_StrangSplitting.hpp>


namespace sweet
{

/*
 * This is a class to register all time steppers
 */
class DESolver_TimeStepperRegistryAll
{
public:
	DESolver_TimeStepperRegistryAll()
	{
	}

	void registerAll(
			sweet::DESolver_TimeTreeNode_Registry &o_timeStepper_registry
	)
	{
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_ExplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_ImplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_AddDETerms>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_StrangSplitting>();

	}
};

}
#endif /* SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_TIMESTEPPERREGISTRY_HPP_ */
