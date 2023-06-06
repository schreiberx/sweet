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
	public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_NegTendencies>
{
public:
	DESolver_TimeStepper_NegTendencies()
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

	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() == 0)
			return error.set("No DE terms specified for time stepper"+getNewLineDebugMessage());

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
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION:
				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				return error.set("Key-value not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
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
				"tendencies"
			);

		_tmpDataContainer.resize(1);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_NegTendencies);
	}

private:
	bool _eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		o_U.op_setZero();

		for (std::size_t i = 0; i < _evalFuns.size(); i++)
		{
			evalTimeStepper(
					i,
					i_U,
					*_tmpDataContainer[0],
					i_simulationTime
				);
			o_U.op_addVector(*_tmpDataContainer[0]);
		}

		o_U.op_mulScalar(-1.0);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "Add(" << std::endl;
		std::cout << newPrefix << "  numDETerms: " << _timeTreeNodes.size() << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
