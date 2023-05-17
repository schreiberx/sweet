/*
 * TimeStepperBase.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPER_BASE_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPER_BASE_HPP_

#include <vector>
#include <string>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include "PDESolver_DataContainer_Base.hpp"
#include "PDESolver_TimeStepping_Tree.hpp"


/*
 * WARNING: Do not include it here. It's included at the end due to circular dependency
 * #include "PDESolver_TimeStepping_Assemblation.hpp"
 */
namespace sweet {
	class PDESolver_TimeStepping_Assemblation;
}


namespace sweet
{

class PDESolver_TimeStepper_Base
{
public:
	sweet::ErrorBase error;

	PDESolver_TimeStepper_Base()
	{
	}

	virtual
	~PDESolver_TimeStepper_Base()
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
	std::shared_ptr<PDESolver_TimeStepper_Base> getNewInstance() = 0;


	/*
	 * Shack registration
	 */
	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) = 0;


	/*
	 * Setup potential internal data structures,
	 * e.g., storage space for RK stages
	 */
	virtual
	bool setupDataContainers(const PDESolver_DataContainer_Base &i_u) = 0;

	virtual
	bool setupFunction(
			std::shared_ptr<sweet::PDESolver_TimeStepping_Tree::Function> &i_function,
			sweet::PDESolver_TimeStepping_Assemblation &i_tsAssemblation
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
	void setTimestepSize(double i_dt) = 0;


	/*
	 * Return the time integration
	 */
	virtual
	void eval_timeIntegration(
			PDESolver_DataContainer_Base &i_u,
			PDESolver_DataContainer_Base &o_u,
			double i_simulation_time
		) = 0;
};

}

// Put at the end here to include it
#include "PDESolver_TimeStepping_Assemblation.hpp"

#endif
