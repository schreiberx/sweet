/*
 * PDETermRegistry.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERTIMESTEPPER_PDETERM_REGISTRY_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERTIMESTEPPER_PDETERM_REGISTRY_HPP_

#include <map>
#include <sweet/core/ErrorBase.hpp>
#include "DESolver_DETerm_Base.hpp"


namespace sweet
{

class DESolver_DETerm_Registry
{
public:
	sweet::ErrorBase error;

private:
	std::map<std::string, DESolver_DETerm_Base*> registry;

public:
	DESolver_DETerm_Registry()
	{
	}

public:
	~DESolver_DETerm_Registry()
	{
		clear();
	}

public:
	template <typename T>
	bool registerPDETerm()
	{
		T *t = new T;
		DESolver_DETerm_Base *b = t;

		// get names of all PDE terms
		std::string ts = b->getImplementedPDETerm();

		auto iter = registry.find(ts);

		if (iter != registry.end())
			return error.set("PDE Term with string "+ts+" already registered");

		registry[ts] = b;
		return true;
	}

	/*
	 * Find a time stepper according to its string
	 */
	bool getPDETermInstance(
			const std::string &i_pdeTermString,
			std::shared_ptr<DESolver_DETerm_Base> &o_pde_term_instance
	)
	{
		auto iter = registry.find(i_pdeTermString);

		if (iter == registry.end())
			return error.set("PDE Term with string "+i_pdeTermString+" not found");

		o_pde_term_instance = iter->second->getNewInstance();
		return true;
	}

	void clear()
	{
	    for (auto it = registry.begin(); it != registry.end(); ++it)
	    	delete it->second;
	}
};

}

#endif
