#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_NEGDETERMS_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_NEGDETERMS_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

/*
 * Negate a PDE term (put a minus sign in front of the tendencies)
 */
class DESolver_TimeStepper_NegTendencies	:
		public sweet::DESolver_TimeTreeNode_Base
{
private:
	double _timestep_size;
	double &dt = _timestep_size;

	// DE term to evaluate
	std::vector<std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>> _timeTreeNode;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _tmpDataContainer;

public:
	DESolver_TimeStepper_NegTendencies()	:
		_timestep_size(-1)
	{
		setEvalAvailable("tendencies");
	}

	~DESolver_TimeStepper_NegTendencies()
	{
		clear();
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("neg");
		retval.push_back("NEG");
		return retval;
	}

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		for (auto &i : _timeTreeNode)
		{
			i->shackRegistration(io_shackDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(*i);
		}

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNode.size() == 0)
			return error.set("No DE terms specified for time stepper"+getNewLineDebugMessage());

		for (auto &i : _timeTreeNode)
		{
			if (!i->isEvalAvailable("tendencies"))
				return error.set("tendencies not available in DE term");
		}

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
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION:
				_timeTreeNode.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNode.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
				break;


			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:

				_timeTreeNode.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNode.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
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


	bool setupConfig(
		const sweet::DESolver_Config_Base &i_deConfig
	) override
	{
		for (auto &i : _timeTreeNode)
		{
			i->setupConfig(i_deConfig);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i);
		}

		_tmpDataContainer.resize(1);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deConfig.getNewDataContainerInstance();

		return true;

	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			delete _tmpDataContainer[i];
		_tmpDataContainer.clear();
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getNewInstance()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_NegTendencies);
	}

	void setTimeStepSize(double i_dt)	override
	{
		_timestep_size = i_dt;

		for (auto &i : _timeTreeNode)
		{
			i->setTimeStepSize(i_dt);
		}
	}

	void eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)	override
	{
		o_U.op_setZero();

		for (auto &i : _timeTreeNode)
		{
			i->eval_tendencies(i_U, *_tmpDataContainer[0], i_simulation_time);
			o_U.op_addVector(*_tmpDataContainer[0]);
		}

		o_U.op_mulScalar(-1.0);
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "Add(" << std::endl;
		std::cout << newPrefix << "  numDETerms: " << _timeTreeNode.size() << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
