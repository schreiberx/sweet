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

#include <sweet/timeTree/DESolver_TimeStepper_AddTendencies.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_NegTendencies.hpp>

#include <sweet/timeTree/DESolver_TimeStepper_ExplicitRungeKutta.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_ImplicitRungeKutta.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_StrangSplitting.hpp>

#include <sweet/timeTree/DESolver_TimeStepper_SubCycling.hpp>

#include <sweet/timeTree/DESolver_TimeStepper_Exponential.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_REXI.hpp>
#include <sweet/timeTree/DESolver_TimeStepper_ETDRK.hpp>


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
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_AddTendencies>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_NegTendencies>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_ExplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_ImplicitRungeKutta>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_StrangSplitting>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_SubCycling>();

		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_Exponential>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_REXI>();
		o_timeStepper_registry.registerTimeTreeNode<sweet::DESolver_TimeStepper_ETDRK>();

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(o_timeStepper_registry);
		return true;
	}
};

}
#endif /* SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_TIMESTEPPERREGISTRY_HPP_ */
