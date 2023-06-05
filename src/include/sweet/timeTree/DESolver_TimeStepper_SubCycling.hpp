#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_SUBCYCLING_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_SUBCYCLING_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeInteriorHelper.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_SubCycling	:
	public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_SubCycling>
{
private:
	// Number of subcycling intervals
	int _subCyclingIntervals;


public:
	DESolver_TimeStepper_SubCycling()
	{
		setEvalAvailable("integration");
	}

	~DESolver_TimeStepper_SubCycling()
	{
		clear();
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("subcyc");
		retval.push_back("subcycling");
		retval.push_back("SUBCYC");
		retval.push_back("SUBCYCLING");
		return retval;
	}

	bool _setupArgumentInternals()
	{
		if (_subCyclingIntervals <= 1)
			return error.set("At least one interval for subcycling required"+getNewLineDebugMessage());

		if (_timeTreeNodes.size() == 0)
			return error.set("Subcycling requires one term/function"+getNewLineDebugMessage());

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

				if (_timeTreeNodes.size() >= 1)
					return error.set("Subcycling only supports a single function/DETerm"+a->getNewLineDebugMessage());
					
				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
				return error.set("Key with functions not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "intervals" || a->key == "i")
				{
					a->getValue(_subCyclingIntervals);
					std::cout << _subCyclingIntervals << std::endl;
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
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

		// provide debug message in case that something goes wrong with the arguments
		setDebugMessage(i_function->getDebugMessage());
		return _setupArgumentInternals();
	}


	bool setupConfigAndGetTimeStepperEval(
		const sweet::DESolver_Config_Base &i_deTermConfig,
		const std::string &i_timeStepperEvalName,
		DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
	) override
	{
		_helperSetupConfigAndGetTimeStepperEval(
				i_deTermConfig,
				i_timeStepperEvalName,
				o_timeStepper,
				"integration"
			);
		
		_tmpDataContainer.resize(1);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();
		
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	inline
	void setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		for (auto &i : _timeTreeNodes)
		{
			i->setTimeStepSize(_timestepSize/_subCyclingIntervals);
		}
	}

	void _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		if (_subCyclingIntervals == 1)
		{
			evalTimeStepper(i_U, o_U, i_simulationTime);
			return;
		}

		int stepsTodo = _subCyclingIntervals;

		// 1st step to use temporary buffer
		evalTimeStepper(i_U, *_tmpDataContainer[0], i_simulationTime);
		stepsTodo--;

		// perform always 2 steps
		for (int i = 1; i < _subCyclingIntervals-1; i+=2)
		{
			evalTimeStepper(*_tmpDataContainer[0], o_U, i_simulationTime);
			evalTimeStepper(o_U, *_tmpDataContainer[0], i_simulationTime);
			stepsTodo -= 2;
		}

		assert(stepsTodo > 0 && stepsTodo <= 2);

		evalTimeStepper(*_tmpDataContainer[0], o_U, i_simulationTime);
		stepsTodo--;

		if (stepsTodo == 0)
			return;

		_tmpDataContainer[0]->swap(o_U);
		evalTimeStepper(*_tmpDataContainer[0], o_U, i_simulationTime);
		stepsTodo--;

		assert(stepsTodo == 0);
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "SUBCYCLING(" << std::endl;
		std::cout << newPrefix << "  subCyclingIntervals: " << _subCyclingIntervals << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
