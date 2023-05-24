/*
 * PDETermRegistry.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPER_REGISTRY_HPP_
#define SRC_PROGRAMS_SIMDATA_PDESOLVER_TIMESTEPPER_REGISTRY_HPP_

#include <map>
#include <sweet/core/ErrorBase.hpp>
#include "DESolver_TimeTreeNode_Base.hpp"


namespace sweet
{

class DESolver_TimeTreeNode_Registry
{
public:
	sweet::ErrorBase error;

private:
	std::map<std::string, std::shared_ptr<DESolver_TimeTreeNode_Base>> registry;

public:
	DESolver_TimeTreeNode_Registry()
	{
	}

public:
	~DESolver_TimeTreeNode_Registry()
	{
		clear();
	}

public:
	template <typename T>
	bool registerTimeTreeNode()
	{
		std::shared_ptr<DESolver_TimeTreeNode_Base> b = std::make_shared<T>();

		// get names of all time steppers
		const std::vector<std::string> ts = b->getNodeNames();

		// insert as keys to registry
		for (auto &iter : ts)
		{
			auto i = registry.find(iter);

			if (i != registry.end())
			{
				return error.set("Time tree node handler for '"+iter+"' already registered!");
			}

			registry[iter] = b;
		}

		return true;
	}

	/*
	 * Find a time stepper according to its string
	 */
	bool getTimeTreeNodeNewInstance(
			const std::string &i_ts_string,
			std::shared_ptr<DESolver_TimeTreeNode_Base> &o_timestepper_instance
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
