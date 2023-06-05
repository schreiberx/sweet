#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_STRANG_SPLITTING_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_STRANG_SPLITTING_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeInteriorHelper.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_StrangSplitting	:
		public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_StrangSplitting>
{
private:
	// Order of Strang splitting
	int _order;

public:
	DESolver_TimeStepper_StrangSplitting()
	{
		setEvalAvailable("integration");
	}

	~DESolver_TimeStepper_StrangSplitting()
	{
		clear();
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("ss");
		retval.push_back("SS");
		return retval;
	}

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		for (auto &i : _timeTreeNodes)
		{
			i->shackRegistration(io_shackDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*i);
		}

		return true;
	}


	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 2)
			return error.set("Only two time steppers supported in this strang splitting"+getNewLineDebugMessage());

		if (_order < 1 || _order > 2)
			return error.set("Only order 1 or 2 allowed for SS method"+getNewLineDebugMessage());

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
				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
				return error.set("Key with functions not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "order" || a->key == "o")
				{
					a->getValue(_order);
					std::cout << _order << std::endl;
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
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

		if (_order == 1)
			_tmpDataContainer.resize(1);
		else if (_order == 2)
			_tmpDataContainer.resize(2);

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		return true;
	}

	inline
	void setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		if (_order == 1)
		{
			for (auto &i : _timeTreeNodes)
			{
				i->setTimeStepSize(i_dt);
			}
			return;
		}

		if (_order == 2)
		{
			_timeTreeNodes[0]->setTimeStepSize(i_dt*0.5);
			_timeTreeNodes[1]->setTimeStepSize(i_dt);
			return;
		}

		SWEETError("Internal error");
	}




	void _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		if (_order == 1)
		{
			_timeTreeNodes[0]->_eval_integration(i_U, *_tmpDataContainer[0], i_simulationTime);
			_timeTreeNodes[1]->_eval_integration(*_tmpDataContainer[0], o_U, i_simulationTime);
			return;
		}
		
		if (_order == 2)
		{
			_timeTreeNodes[0]->_eval_integration(i_U, *_tmpDataContainer[0], i_simulationTime);
			_timeTreeNodes[1]->_eval_integration(*_tmpDataContainer[0], *_tmpDataContainer[1], i_simulationTime);
			_timeTreeNodes[0]->_eval_integration(*_tmpDataContainer[1], o_U, i_simulationTime+0.5*dt);
			return;
		}
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "SS(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
