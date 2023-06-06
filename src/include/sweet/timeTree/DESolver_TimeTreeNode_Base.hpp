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
	 * Type define 'EvalFun' which is required to time step
	 */
	typedef bool (DESolver_TimeTreeNode_Base::*EvalFun)(
					const DESolver_DataContainer_Base &i_U,
					DESolver_DataContainer_Base &o_U,
					double i_simulationTime
			);

	/*
	 * Keep track of a list of implemented evaluation functions
	 */
protected:
	std::vector<std::string> _registeredEvalFunctions;

public:
	DESolver_TimeTreeNode_Base()
	{
	}

	virtual
	~DESolver_TimeTreeNode_Base()
	{
	}

	/**
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


	/**
	 * Return a new instance of this class
	 */
	virtual
	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew() = 0;


	/**
	 * Return a copy of this class
	 *
	 * This is a special case where we might like to reuse the already utilized shacks
	 */
	virtual
	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy() = 0;


	/**
	 * Shack registration
	 */
	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) = 0;


	/**
	 * Setup the treenode function
	 */
	virtual
	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::DESolver_TimeStepping_Tree::Function> &i_function,
			sweet::DESolver_TimeStepping_Assemblation &i_tsAssemblation
	)
	{
		return error.set("setupTreeNodeByFunction not supported for this node");
	}


	/**
	 * Setup operators and potential internal data structures,
	 * e.g., storage space for RK stages
	 */
	virtual
	bool setupConfigAndGetTimeStepperEval(
			const sweet::DESolver_Config_Base &i_deTermConfig,
			const std::string &i_timeStepperEvalName,
			EvalFun &o_timeStepper
		) = 0;


	/**
	 * This allows customizing the time stepper by its parent class
	 *
	 * This can be helpful, e.g., for exponential integration to specify the particular phi function:
	 *
	 * "expIntegratorFunction" => "phi0"
	 * or
	 * "expIntegratorFunction" => "phi1"
	 * ...
	 */
	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	){
		return error.set("setupByKeyValue() not available in tree node '"+getNodeNames()[0]+"'");
	};

public:
	bool _helperGetTimeStepperEval(
			const std::string &i_timeStepperEvalName,
			EvalFun &o_timeStepper
	)
	{
		if (!isEvalAvailable(i_timeStepperEvalName))
			return error.set("Time stepper evaluation '"+i_timeStepperEvalName+"' not available or not registered");

		if (i_timeStepperEvalName == "integration")
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_integration;
			return true;
		}

		if (i_timeStepperEvalName == "tendencies")
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_tendencies;
			return true;
		}

		if (i_timeStepperEvalName == "eulerBackward")
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_eulerBackward;
			return true;
		}

		if (i_timeStepperEvalName == "eulerForward")
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_eulerForward;
			return true;
		}

		if (i_timeStepperEvalName == "exponential")
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_exponential;
			return true;
		}

		return error.set("Time stepper '"+i_timeStepperEvalName+"' not found");
	}


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
private:
	bool isEvalAvailable(const std::string &i_evalFunction)
	{
		for (auto &e: _registeredEvalFunctions)
		{
			if (e == i_evalFunction)
				return true;
		}

		return false;
	}


	/**
	 * Set an evaluation function to be available
	 */
public:
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

		_registeredEvalFunctions.push_back(i_evalFunction);
	}


	/*
	 * Return the time integration
	 */
	virtual
	bool _eval_integration(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	) {
		return error.set("_eval_integration() not available");
	};


	/*
	 * Optional: Return the time tendencies of the DE term
	 */
public:
	virtual
	bool _eval_tendencies(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
	){
		return error.set("_eval_tendencies() not available");
	};

	/*
	 * Optional: Return the forward Euler evaluation of the term:
	 *
	 * (U^{n+1}-U^{n})/Dt = dU^{n}/dt
	 * <=> U^{n+1} = U^n + Dt* (dU^{n}n/dt)
	 */
public:
	virtual
	bool _eval_eulerForward(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
	){
		return error.set("_eval_eulerForward() not available");
	};

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
	bool _eval_eulerBackward(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
	){
		return error.set("_eval_eulerBackward() not available");
	};

	/*
	 * Optional: Return an evaluation of the exponential term
	 *
	 * U^{n+1} = exp(Dt*L) U^n
	 */
	virtual
	bool _eval_exponential(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
	){
		return error.set("_eval_exponential() not available");
	};

	/*!
	 * Message to help debugging errors in the time tree (e.g. about the parsing of the time tree string)
	 */
	std::string _debugMessage;

	/*!
	 * Set the debug message
	 */
	std::string setDebugMessage(
			const std::string &i_debugMessage	//!< Message for debugging
	)
	{
		_debugMessage = i_debugMessage;
		return _debugMessage;
	}

	/*!
	 * Return the debug message including a newline.
	 * This is handy for figuring out wrong parameters in the time tree.
	 *
	 * \return String with hopefully helpful debug information
	 */
	std::string getNewLineDebugMessage()
	{
		if (_debugMessage != "")
			return "\n"+_debugMessage;

		return _debugMessage;
	}
};

}

#endif
