/*
 * TimeStepperBase.hpp
 */

#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERMBASE_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERMBASE_HPP_

#include <vector>
#include <string>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/ErrorBase.hpp>
#include "DESolver_DETerm_Base.hpp"
#include "DESolver_DataContainer_Base.hpp"
#include "DESolver_Config_Base.hpp"


namespace sweet
{

class DESolver_DETerm_Base
{
public:
	ErrorBase error;

	DESolver_DETerm_Base()
	{
	}

	virtual
	~DESolver_DETerm_Base()
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
	std::shared_ptr<DESolver_DETerm_Base> getNewInstance() = 0;

	/*
	 * Shack registration
	 */
	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) = 0;


#if 0
	/*
	 * Setup potential internal data structures
	 */
	virtual
	bool setupOpsAndDataContainers(
			const sweet::SphereOperators *io_ops,
			const DESolver_DataContainer_Base &i_u
		) = 0;
#endif

	virtual
	bool setupDETermConfig(
		const sweet::DESolver_Config_Base &i_deTermConfig
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
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
		){};

	/*
	 * Optional: Return the forward Euler evaluation of the term:
	 *
	 * (U^{n+1}-U^{n})/Dt = dU^{n}/dt
	 * <=> U^{n+1} = U^n + Dt* (dU^{n}n/dt)
	 */
	virtual
	void eval_eulerForward(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
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
	void eval_eulerBackward(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
		){};

	/*
	 * Optional: Return an evaluation of the exponential term
	 *
	 * U^{n+1} = exp(Dt*L) U^n
	 */
	virtual
	void eval_exponential(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
		){};


	/**
	 * Check whether some evaluation is available.
	 *
	 * This only works with a hack in the GCC compiler.
	 *
	 * For other compilers it will always return that the function is available.
	 */

	bool isEvalAvailable(
			const std::string &i_functionName
	)
	{
#ifndef __GNUC__
		return true;
#else
	#ifdef __clang__
		return true;
	#else

		#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wpmf-conversions"


		if (i_functionName == "eval_tendencies")
		{
			void *a = (void*)(this->*(&DESolver_DETerm_Base::eval_tendencies));
			void *b = (void*)&DESolver_DETerm_Base::eval_tendencies;
			return a != b;
		}
		if (i_functionName == "eval_eulerForward")
		{
			void *a = (void*)(this->*(&DESolver_DETerm_Base::eval_eulerForward));
			void *b = (void*)&DESolver_DETerm_Base::eval_eulerForward;
			return a != b;
		}
		if (i_functionName == "eval_eulerBackward")
		{
			void *a = (void*)(this->*(&DESolver_DETerm_Base::eval_eulerBackward));
			void *b = (void*)&DESolver_DETerm_Base::eval_eulerBackward;
			return a != b;
		}
		if (i_functionName == "eval_exponential")
		{
			void *a = (void*)(this->*(&DESolver_DETerm_Base::eval_exponential));
			void *b = (void*)&DESolver_DETerm_Base::eval_exponential;
			return a != b;
		}

		return error.set("Invalid function '"+i_functionName+"'");

		#pragma GCC diagnostic pop

		return true;
	#endif
#endif
	}
	bool isEvalAvailable(
			const DESolver_DETerm_Base &i_deTermBase,
			const std::string &i_functionName
	)
	{
#ifndef __GNUC__
		return true;
#else
	#ifdef __clang__
		return true;
	#else

		#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wpmf-conversions"


		if (i_functionName == "eval_tendencies")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_DETerm_Base::eval_tendencies));
			void *b = (void*)&DESolver_DETerm_Base::eval_tendencies;
			return a != b;
		}
		if (i_functionName == "eval_eulerForward")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_DETerm_Base::eval_eulerForward));
			void *b = (void*)&DESolver_DETerm_Base::eval_eulerForward;
			return a != b;
		}
		if (i_functionName == "eval_eulerBackward")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_DETerm_Base::eval_eulerBackward));
			void *b = (void*)&DESolver_DETerm_Base::eval_eulerBackward;
			return a != b;
		}
		if (i_functionName == "eval_exponential")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_DETerm_Base::eval_exponential));
			void *b = (void*)&DESolver_DETerm_Base::eval_exponential;
			return a != b;
		}

		return error.set("Invalid function '"+i_functionName+"'");

		#pragma GCC diagnostic pop

		return true;
	#endif
#endif
	}
};

}

#endif
