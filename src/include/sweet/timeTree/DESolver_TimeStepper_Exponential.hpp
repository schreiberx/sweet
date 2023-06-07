#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPONENTIAL_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPONENTIAL_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeInteriorHelper.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_Exponential	:
	public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_Exponential>
{
private:
	std::string _expFunctionString;

public:
	DESolver_TimeStepper_Exponential()
	{
		setEvalAvailable(EVAL_EXPONENTIAL);
		setEvalAvailable(EVAL_INTEGRATION);
	}


	~DESolver_TimeStepper_Exponential()
	{
		clear();
	}

	DESolver_TimeStepper_Exponential(
			const DESolver_TimeStepper_Exponential &i_src
	)	:
		DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_Exponential>(i_src)
	{
		_expFunctionString = i_src._expFunctionString;
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("exp");
		retval.push_back("EXP");
		return retval;
	}


	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() == 0)
			return error.set("Some time node term needs to be given"+getNewLineDebugMessage());

		return true;
	}


	virtual
	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::DESolver_TimeStepping_Tree::Function> &i_function,
			sweet::DESolver_TimeStepping_Assemblation &i_tsAssemblation
	)	override
	{
		for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
		{
			sweet::DESolver_TimeStepping_Tree::Argument *a = iter->get();

			switch(a->argType)
			{
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION:
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
				if (_timeTreeNodes.size() == 1)
					return error.set("Only one term for exponential integration supported");

				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);

				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() != 0)
					return error.set("Only one term for exponential integration supported");

				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			default:
				SWEETError("Internal error");
				return error.set("Internal error");
			}
		}

		// Provide debug message in case that something goes wrong with the arguments
		setDebugMessage(i_function->getDebugMessage());
		return _setupArgumentInternals();
	}

	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override
	{
		if (i_key == "expIntegrationFunction")
		{
			assert(_timeTreeNodes[0] != nullptr);
			_timeTreeNodes[0]->setupByKeyValue(i_key, i_value);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

			_expFunctionString = i_value;
			return true;
		}

		return false;
	}

	bool setupConfigAndGetTimeStepperEval(
		const sweet::DESolver_Config_Base &i_deTermConfig,
		EVAL_TYPES i_evalType,
		DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
	) override
	{
		_helperSetupConfigAndGetTimeStepperEval(
				i_deTermConfig,
				i_evalType,
				o_timeStepper,
				EVAL_EXPONENTIAL
			);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}


#if 0
	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_Exponential);
	}
#endif

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_Exponential(*this));
	}


private:
	bool _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		assert(_timeTreeNodes[0] != nullptr);
		return evalTimeStepper(i_U, o_U, i_simulationTime);
	}



	/*
	 * We also provide an exponential time integration for this one
	 * in order to transparently support 'exponential' time integration for
	 * either DE terms themselves, EXP and also REXI evaluations.
	 */
private:
	bool _eval_exponential(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		assert(_timeTreeNodes[0] != nullptr);
		return evalTimeStepper(i_U, o_U, i_simulationTime);
	}


	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "EXP(" << std::endl;
		std::cout << newPrefix << "  expFunctionString: '" << _expFunctionString << "'" << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
