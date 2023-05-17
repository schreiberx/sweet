/*
 * PDETermRegistry.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPER_REGISTRY_HPP_
#define SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPER_REGISTRY_HPP_

#include <map>
#include <sweet/core/ErrorBase.hpp>
#include "PDESolver_TimeStepper_Base.hpp"


namespace sweet
{

class PDESolver_TimeStepper_Registry
{
public:
	sweet::ErrorBase error;

private:
	std::map<std::string, std::shared_ptr<PDESolver_TimeStepper_Base>> registry;

public:
	PDESolver_TimeStepper_Registry()
	{
	}

public:
	~PDESolver_TimeStepper_Registry()
	{
		clear();
	}

public:
	template <typename T>
	bool registerTimeStepper()
	{
		std::shared_ptr<PDESolver_TimeStepper_Base> b = std::make_shared<T>();

		// get names of all time steppers
		const std::vector<std::string> ts = b->getImplementedTimeSteppers();

		// insert as keys to registry
		for (std::size_t i = 0; i < ts.size(); i++)
		{
			auto iter = registry.find(ts[i]);

			if (iter != registry.end())
				return error.set("Time stepper with string "+ts[i]+" already registered");

			registry[ts[i]] = b;
		}

		return true;
	}

	/*
	 * Find a time stepper according to its string
	 */
	bool getTimeStepperNewInstance(
			const std::string &i_ts_string,
			std::shared_ptr<PDESolver_TimeStepper_Base> &o_timestepper_instance
	)
	{
		auto iter = registry.find(i_ts_string);

		if (iter == registry.end())
			return error.set("Timestepper with string "+i_ts_string+" not found");

		o_timestepper_instance = iter->second->getNewInstance();
		return true;
	}

	void clear()
	{
		registry.clear();
	}
};

}

#endif
