/*
 * TimeStepperBase.hpp
 */

/*
 * Make sure that we always try to include these classes to due
 * forward declarations
 */
#include "DESolver_TimeStepping_Assemblation.hpp"

#ifndef SRC_SWEET_TIMESTEPPER_BASE_HPP_
#define SRC_SWEET_TIMESTEPPER_BASE_HPP_

#include <vector>
#include <string>
#include <complex>
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


	/*!
	 * Different evaluations for time integration
	 */
	enum EVAL_TYPES
	{
		EVAL_INTEGRATION,
		EVAL_TENDENCIES,
		EVAL_EULER_BACKWARD,
		EVAL_REXI_TERM,
		EVAL_EULER_FORWARD,
		EVAL_EXPONENTIAL
	};


	/*!
	 * Type define 'EvalFun' which is required to time step
	 */
	typedef bool (DESolver_TimeTreeNode_Base::*EvalFun)(
					const DESolver_DataContainer_Base &i_U,
					DESolver_DataContainer_Base &o_U,
					double i_simulationTime
			);

	/*!
	 * Keep track of a list of implemented evaluation functions
	 */
protected:
	std::vector<EVAL_TYPES> _registeredEvalTypes;

public:
	DESolver_TimeTreeNode_Base()
	{
	}

	virtual
	~DESolver_TimeTreeNode_Base()
	{
	}

	/*!
	 * Copy constructor
	 */
public:
	DESolver_TimeTreeNode_Base(
			const DESolver_TimeTreeNode_Base &i_src
	)
	{
		_registeredEvalTypes = i_src._registeredEvalTypes;
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

#if 0
	/**
	 * Return a new instance of this class
	 */
	virtual
	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew() = 0;
#endif

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
			EVAL_TYPES i_evalType,
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
		return error.set("setupByKeyValue(std::string,std::string) not available in tree node '"+getNodeNames()[0]+"'");
	};


	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
	){
		return error.set("setupByKeyValue(std::string,std::complex<double>) not available in tree node '"+getNodeNames()[0]+"'");
	};

public:
	bool _helperGetTimeStepperEval(
			EVAL_TYPES i_evalType,
			EvalFun &o_timeStepper
	)
	{
		if (!isEvalAvailable(i_evalType))
		{
			std::ostringstream oss;
			oss << "Time stepper evaluation '" << i_evalType << "' not available or not registered";
			return error.set(oss.str());
		}

		if (i_evalType == EVAL_INTEGRATION)
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_integration;
			return true;
		}

		if (i_evalType == EVAL_TENDENCIES)
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_tendencies;
			return true;
		}

		if (i_evalType == EVAL_EULER_BACKWARD)
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_eulerBackward;
			return true;
		}

		if (i_evalType == EVAL_REXI_TERM)
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_rexiTerm;
			return true;
		}

		if (i_evalType == EVAL_EULER_FORWARD)
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_eulerForward;
			return true;
		}

		if (i_evalType == EVAL_EXPONENTIAL)
		{
			o_timeStepper = &DESolver_TimeTreeNode_Base::_eval_exponential;
			return true;
		}

		std::ostringstream oss;
		oss << "Invalid Time stepper evaluation '" << i_evalType << "'";
		return error.set(oss.str());
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
	bool isEvalAvailable(EVAL_TYPES &i_evalType)
	{
		for (auto &e: _registeredEvalTypes)
		{
			if (e == i_evalType)
				return true;
		}

		return false;
	}


	/**
	 * Set an evaluation function to be available
	 */
public:
	void setEvalAvailable(EVAL_TYPES i_evalType)
	{
		_registeredEvalTypes.push_back(i_evalType);
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
	 * Optional: Return the backward Euler evaluation with complex fields
	 */
	virtual
	bool _eval_rexiTerm(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_time_stamp
	){
		return error.set("_eval_rexiTerm() not available");
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
