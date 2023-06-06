#ifndef SRC_PROGRAMS_SIMDATA_PDESWESPHERE_DESOLVER_CONFIG_HPP_
#define SRC_PROGRAMS_SIMDATA_PDESWESPHERE_DESOLVER_CONFIG_HPP_

#include <sweet/timeTree/DESolver_Config_Base.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperatorsComplex.hpp>
#include "PDESWESphere_DataContainer.hpp"

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
	sweet::SphereOperators *ops;
	sweet::SphereOperatorsComplex *opsComplex;

	/*
	 * Return a new instance of a data container.
	 *
	 * This is what will be used by time steppers and/or DE term implementations
	 */
	sweet::DESolver_DataContainer_Base* getNewDataContainerInstance(
			int i_id = -1
	) const override
	{
		PDESWESphere_DataContainer *retval = new PDESWESphere_DataContainer;
		retval->setup_like(*myDataContainer);

		return retval;
	}
};



#endif /* SRC_PROGRAMS_SIMDATA_PDESWESPHERE_DESOLVER_CONFIG_HPP_ */
