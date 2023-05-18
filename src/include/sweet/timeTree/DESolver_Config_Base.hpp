#ifndef SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_DETERMCONFIG_BASE_HPP_
#define SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_DETERMCONFIG_BASE_HPP_


#include "DESolver_DataContainer_Base.hpp"


namespace sweet
{

class DESolver_Config_Base
{
public:
	DESolver_Config_Base()
	{
	}

	virtual ~DESolver_Config_Base()
	{
	}

	/**
	 * Return a new instance of a data object related to this time stepper
	 */
	virtual
	DESolver_DataContainer_Base* getNewDataContainerInstance(int i_id = -1) const = 0;
};

}

#endif /* SRC_INCLUDE_SWEET_TIMENEW_DESOLVER_DETERMCONFIG_BASE_HPP_ */
