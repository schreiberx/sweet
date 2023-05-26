#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPONENTIAL_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPONENTIAL_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_Exponential	:
		public sweet::DESolver_TimeTreeNode_Base
{
private:
	double _timestep_size;
	double &dt = _timestep_size;

	// Order of Strang splitting
	int _order;

	// DE term to evaluate
	std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> _timeTreeNode;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _tmpDataContainer;


public:
	DESolver_TimeStepper_Exponential()	:
		_timestep_size(-1)
	{
		setEvalAvailable("integration");
	}


	~DESolver_TimeStepper_Exponential()
	{
		clear();
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("exp");
		retval.push_back("EXP");
		return retval;
	}


	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		if (_timeTreeNode != nullptr)
		{
			_timeTreeNode->shackRegistration(io_shackDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(*_timeTreeNode);
		}

		return true;
	}


	bool _setupArgumentInternals()
	{
		if (_timeTreeNode == nullptr)
			return error.set("Some time node term needs to be given"+getNewLineDebugMessage());

		if (!_timeTreeNode->isEvalAvailable("exponential"))
			return error.set("exponential evaluation not available in DE term for exponential integrator");

		return true;
	}


	virtual
	bool setupFunction(
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
				if (_timeTreeNode != nullptr)
					return error.set("Only one term for exponential integration supported");

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNode
					);

				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
				return error.set("Key with functions not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNode != nullptr)
					return error.set("Only one term for exponential integration supported");

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNode
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


	bool setupConfig(
		const sweet::DESolver_Config_Base &i_deConfig
	) override
	{
		_timeTreeNode->setupConfig(i_deConfig);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNode);

		_tmpDataContainer.resize(0);

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deConfig.getNewDataContainerInstance();

		if (!_timeTreeNode->isEvalAvailable("exponential"))
			return error.set("exponential integration not available in exp term");

		return true;

	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			delete _tmpDataContainer[i];

		_tmpDataContainer.clear();

		_timeTreeNode.reset();
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getNewInstance()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_Exponential);
	}

	void setTimeStepSize(double i_dt)	override
	{
		_timestep_size = i_dt;

		_timeTreeNode->setTimeStepSize(_timestep_size);
	}

	void eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)	override
	{
		assert(_timeTreeNode != nullptr);
		_timeTreeNode->eval_exponential(i_U, o_U, i_simulation_time);
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "EXP(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
