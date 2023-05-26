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
#include "DESolver_Config_Base.hpp"


/*
 * WARNING: Do not include it here. It's included at the end due to circular dependency
 * #include "DESolver_TimeStepping_Assemblation.hpp"
 */
namespace sweet {
	class DESolver_TimeStepping_Assemblation;
}


namespace sweet
{

class DESolver_TimeTreeNode_Base
{
public:
	sweet::ErrorBase error;

	/*
	 * Keep track of a list of implemented evaluation functions
	 */
private:
	std::vector<std::string> existingEvalFunctions;

public:
	DESolver_TimeTreeNode_Base()
	{
	}

	virtual
	~DESolver_TimeTreeNode_Base()
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
	getNodeNames() = 0;

	virtual
	std::shared_ptr<DESolver_TimeTreeNode_Base> getNewInstance() = 0;


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
	)
	{
		return error.set("setupFunction not supported for this node");
	}


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
	 * Return true if an evaluation function exists.
	 */
	bool isEvalAvailable(const std::string &i_evalFunction)
	{
		for (auto &e: existingEvalFunctions)
		{
			if (e == i_evalFunction)
				return true;
		}

		return false;
	}


	/*
	 * Set an evaluation function to be available
	 */
	void setEvalAvailable(const std::string &i_evalFunction)
	{
		if (
			i_evalFunction != "integration" &&
			i_evalFunction != "tendencies" &&
			i_evalFunction != "eulerForward" &&
			i_evalFunction != "eulerBackward" &&
			i_evalFunction != "exponential"
		)
			SWEETError("Evaluation function '"+i_evalFunction+"' doesn't exist");

		existingEvalFunctions.push_back(i_evalFunction);
	}


	/*
	 * Return the time integration
	 */
	virtual
	void eval_integration(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulation_time
		) {};


	/*
	 * Optional: Return the time tendencies of the DE term
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

#if 0
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


		if (i_functionName == "tendencies")
		{
			void *a = (void*)(this->*(&DESolver_TimeTreeNode_Base::eval_tendencies));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_tendencies;
			return a != b;
		}
		if (i_functionName == "eulerForward")
		{
			void *a = (void*)(this->*(&DESolver_TimeTreeNode_Base::eval_eulerForward));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_eulerForward;
			return a != b;
		}
		if (i_functionName == "eulerBackward")
		{
			void *a = (void*)(this->*(&DESolver_TimeTreeNode_Base::eval_eulerBackward));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_eulerBackward;
			return a != b;
		}
		if (i_functionName == "exponential")
		{
			void *a = (void*)(this->*(&DESolver_TimeTreeNode_Base::eval_exponential));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_exponential;
			return a != b;
		}
		if (i_functionName == "timeIntegration")
		{
			void *a = (void*)(this->*(&DESolver_TimeTreeNode_Base::eval_integration));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_integration;
			return a != b;
		}

		return error.set("Invalid function '"+i_functionName+"'");

		#pragma GCC diagnostic pop

		return true;
	#endif
#endif
	}

	bool isEvalAvailable(
			const DESolver_TimeTreeNode_Base &i_deTermBase,
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


		if (i_functionName == "tendencies")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_TimeTreeNode_Base::eval_tendencies));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_tendencies;
			return a != b;
		}
		if (i_functionName == "eulerForward")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_TimeTreeNode_Base::eval_eulerForward));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_eulerForward;
			return a != b;
		}
		if (i_functionName == "eulerBackward")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_TimeTreeNode_Base::eval_eulerBackward));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_eulerBackward;
			return a != b;
		}
		if (i_functionName == "exponential")
		{
			void *a = (void*)(i_deTermBase.*(&DESolver_TimeTreeNode_Base::eval_exponential));
			void *b = (void*)&DESolver_TimeTreeNode_Base::eval_exponential;
			return a != b;
		}

		return error.set("Invalid function '"+i_functionName+"'");

		#pragma GCC diagnostic pop

		return true;
	#endif
#endif
	}
#endif

	std::string _debugMessage;

	std::string setDebugMessage(const std::string &i_debugMessage)
	{
		_debugMessage = i_debugMessage;
		return _debugMessage;
	}
	std::string getNewLineDebugMessage()
	{
		if (_debugMessage != "")
			return "\n"+_debugMessage;

		return _debugMessage;
	}
};

}

#endif
