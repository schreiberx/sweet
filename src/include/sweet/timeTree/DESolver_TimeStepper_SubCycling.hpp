#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_SUBCYCLING_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_SUBCYCLING_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_SubCycling	:
		public sweet::DESolver_TimeTreeNode_Base
{
private:
	double _timestep_size;
	double &dt = _timestep_size;

	// Number of subcycling intervals
	int _subCyclingIntervals;

	// DE term to evaluate
	std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> _timeTreeNode;
	DESolver_TimeTreeNode_Base::EvalFun _evalFun;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _tmpDataContainer;

public:
	DESolver_TimeStepper_SubCycling()	:
		_timestep_size(-1)
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

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		_timeTreeNode->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(*_timeTreeNode);

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (_subCyclingIntervals <= 1)
			return error.set("At least one interval for subcycling required"+getNewLineDebugMessage());

		if (_timeTreeNode == nullptr)
			return error.set("Subcycling requires one term/function"+getNewLineDebugMessage());

		for (auto &i : _timeTreeNode)
		{
			if (!i->isEvalAvailable("integration"))
				return error.set("eval of 'integration' missing in term '"+i->getNodeNames()[0]+"'");
		}

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

				if (_timeTreeNode.size() >= 1)
					return error.set("Subcycling only supports a single function/DETerm"+a->getNewLineDebugMessage());

				_timeTreeNode.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNode.back()
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
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				_timeTreeNode.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNode.back()
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
		for (auto &i : _timeTreeNode)
		{
			i->setupConfig(i_deTermConfig);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i);
		}

		if (!_timeTreeNode[0]->isEvalAvailable("integration"))
			return error.set("integration not available in 1st SS term");

		// default setup
		DESolver_TimeTreeNode_Base::_helperSetupConfigAndGetTimeStepperEval(
				i_timeStepperEvalName,
				o_timeStepper
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;

	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			delete _tmpDataContainer[i];
		_tmpDataContainer.clear();

		_timeTreeNode.clear();
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getNewInstance()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_SubCycling);
	}

	void setTimeStepSize(double i_dt)	override
	{
		_timestep_size = i_dt;

		_timeTreeNode[0]->setTimeStepSize(dt/_subCyclingIntervals);
	}

	void _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)	override
	{
		if (_subCyclingIntervals == 1)
		{
			_timeTreeNode[0]->_eval_integration(i_U, *_tmpDataContainer[0], i_simulation_time);
			_timeTreeNode[1]->_eval_integration(*_tmpDataContainer[0], o_U, i_simulation_time);
			return;
		}
		
		if (_subCyclingIntervals == 2)
		{
			_timeTreeNode[0]->_eval_integration(i_U, *_tmpDataContainer[0], i_simulation_time);
			_timeTreeNode[1]->_eval_integration(*_tmpDataContainer[0], *_tmpDataContainer[1], i_simulation_time);
			_timeTreeNode[0]->_eval_integration(*_tmpDataContainer[1], o_U, i_simulation_time+0.5*dt);
			return;
		}
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "SS(" << std::endl;
		std::cout << newPrefix << "  order: " << _subCyclingIntervals << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
