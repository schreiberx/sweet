/*
 * TimeStepperBase.hpp
 */

/*
 * Make sure that we always try to include these classes to due
 * forward declarations
 */
#include "DESolver_TimeStepping_Assemblation.hpp"

#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPER_BASE_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPER_BASE_HPP_

#include <vector>
#include <string>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include "DESolver_DataContainer_Base.hpp"
#include "DESolver_TimeStepping_Tree.hpp"


/*
 * WARNING: Do not include it here. It's included at the end due to circular dependency
 * #include "DESolver_TimeStepping_Assemblation.hpp"
 */
namespace sweet {
	class DESolver_TimeStepping_Assemblation;
}


namespace sweet
{

class DESolver_TimeStepper_Base
{
public:
	sweet::ErrorBase error;

	DESolver_TimeStepper_Base()
	{
	}

	virtual
	~DESolver_TimeStepper_Base()
	{
	}

	/*
	 * Return string of implemented PDE term.
	 *
	 * e.g.
	 * 	'ERK': fast gravity modes
	 * 	'SS': coriolis effect
	 * 	'l': lg+lc
	 */
	virtual
	const std::vector<std::string>
	getImplementedTimeSteppers() = 0;

	virtual
	std::shared_ptr<DESolver_TimeStepper_Base> getNewInstance() = 0;


	/*
	 * Shack registration
	 */
	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) = 0;


	virtual
	bool setupFunction(
			std::shared_ptr<sweet::DESolver_TimeStepping_Tree::Function> &i_function,
			sweet::DESolver_TimeStepping_Assemblation &i_tsAssemblation
	) = 0;


#if 0
	/*
	 * Setup operators and potential internal data structures,
	 * e.g., storage space for RK stages
	 */
	virtual
	bool setupOpsAndDataContainers(
		const sweet::SphereOperators *io_ops,
		const sweet::DESolver_DataContainer_Base &i_U
	) = 0;
#endif


	/*
	 * Setup operators and potential internal data structures,
	 * e.g., storage space for RK stages
	 */
	virtual
	bool setupConfig(
		const sweet::DESolver_Config_Base &i_deTermConfig
	) = 0;

	/*
	 * Cleanup internal data structures
	 */
	virtual
	void clear() = 0;

	/*
	 * Set the time step size Dt.
	 *
	 * This is required, e.g., to setup certain data structures for an implicit time steppers
	 */
	virtual
	void setTimeStepSize(double i_dt) = 0;


	/*
	 * Return the time integration
	 */
	virtual
	void eval_timeIntegration(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulation_time
		) = 0;
};

}

#endif
