/*
 * TimeStepperBase.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERMBASE_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERMBASE_HPP_

#include <vector>
#include <string>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/ErrorBase.hpp>
#include "PDESolver_PDETerm_Base.hpp"
#include "PDESolver_DataContainer_Base.hpp"


namespace sweet
{

class PDESolver_PDETerm_Base
{
public:
	ErrorBase error;

	PDESolver_PDETerm_Base()
	{
	}

	virtual
	~PDESolver_PDETerm_Base()
	{
	}

	/*
	 * Return string of implemented PDE term.
	 *
	 * e.g.
	 * 	'lg': fast gravity modes
	 * 	'lc': coriolis effect
	 * 	'l': lg+lc
	 */
	virtual
	const char* getImplementedPDETerm() = 0;


	/*
	 * Return a new instance
	 */
	virtual
	std::shared_ptr<PDESolver_PDETerm_Base> getNewInstance() = 0;

	/*
	 * Shack registration
	 */
	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) = 0;


	/*
	 * Setup potential internal data structures
	 */
	virtual
	bool setupOpsAndDataContainers(
			const sweet::SphereOperators *io_ops,
			const PDESolver_DataContainer_Base &i_u
		) = 0;


	/*
	 * Set the time step size Dt.
	 *
	 * This is required, e.g., to setup certain data structures for an implicit time steppers
	 */
	virtual
	void setTimestepSize(double i_dt) = 0;

	/*
	 * Optional: Return the time tendencies of the PDE term
	 */
	virtual
	void eval_tendencies(
			const PDESolver_DataContainer_Base &i_U,
			PDESolver_DataContainer_Base &o_U,
			double i_time_stamp
		){};

	/*
	 * Optional: Return the forward Euler evaluation of the term:
	 *
	 * (U^{n+1}-U^{n})/Dt = dU^{n}/dt
	 * <=> U^{n+1} = U^n + Dt* (dU^{n}n/dt)
	 */
	virtual
	void eval_euler_forward(
			const PDESolver_DataContainer_Base &i_U,
			PDESolver_DataContainer_Base &o_U,
			double i_time_stamp
		){};

	/*
	 * Optional: Return the backward Euler evaluation of the term:
	 *
	 * ( U^{n+1} - U^{n} ) / Dt = d/dt U^{n+1}
	 * <=> U^{n+1} - U^{n} = Dt * d/dt U^{n+1}
	 * <=> (I - Dt * d/dt) U^{n+1} = U^{n}
	 *
	 * For a linear operator U_t = LU we would obtain
	 * <=> (I - Dt * L) U^{n+1} = U^{n}
	 */
	virtual
	void eval_euler_backward(
			const PDESolver_DataContainer_Base &i_U,
			PDESolver_DataContainer_Base &o_U,
			double i_time_stamp
		){};

	/*
	 * Optional: Return an evaluation of the exponential term
	 *
	 * U^{n+1} = exp(Dt*L) U^n
	 */
	virtual
	void eval_exponential(
			const PDESolver_DataContainer_Base &i_U,
			PDESolver_DataContainer_Base &o_U,
			double i_time_stamp
		){};
};

}

#endif
