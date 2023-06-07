#ifndef SRC_PROGRAMS_SIMDATA_PDESWESPHERE_DESOLVER_CONFIG_HPP_
#define SRC_PROGRAMS_SIMDATA_PDESWESPHERE_DESOLVER_CONFIG_HPP_

#include <sweet/timeTree/DESolver_Config_Base.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperatorsComplex.hpp>
#include "PDESWESphere_DataContainer.hpp"
//#include "PDESWESphere_DataContainerComplex.hpp"

/*
 * A special class which is forwarded to all
 *
 *  - time stepper instances and
 *  - DE term instances
 *
 * to set up data buffers and other things which are required.
 */
class PDESWESphere_DESolver_Config:
		public sweet::DESolver_Config_Base
{
public:
	/*
	 * Just a pointer to an existing data container
	 */
	const PDESWESphere_DataContainer *myDataContainer;
	//const PDESWESphere_DataContainerComplex *myDataContainerComplex;

	sweet::SphereOperators *ops;
	sweet::SphereOperatorsComplex *opsComplex;


	PDESWESphere_DESolver_Config()
	{
		myDataContainer = nullptr;
		//myDataContainerComplex = nullptr;
		ops = nullptr;
		opsComplex = nullptr;
	}


	/*
	 * Return a new instance of a data container.
	 *
	 * This is what will be used by time steppers and/or DE term implementations
	 */
	sweet::DESolver_DataContainer_Base* getNewDataContainerInstance(
			int i_id = -1
	) const override
	{
		if (i_id == -1)
		{
			/*
			 * Real valued
			 */
			PDESWESphere_DataContainer *retval = new PDESWESphere_DataContainer;
			retval->setup(ops->sphereDataConfig);
			return retval;
		}
#if 0
		if (i_id == 1)
		{
			/*
			 * Complex valued
			 */
			PDESWESphere_DataContainerComplex *retval = new PDESWESphere_DataContainerComplex;
			retval->setup(opsComplex->sphereDataConfig);
			return retval;
		}
#endif
		SWEETError("Invalid data type id");
		return 0;
	}
};



#endif /* SRC_PROGRAMS_SIMDATA_PDESWESPHERE_DESOLVER_CONFIG_HPP_ */
