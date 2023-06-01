#ifndef SRC_SWEET_TIMETREE_TIMESTEPPER_NODE_INTERIOR_HPP_
#define SRC_SWEET_TIMETREE_TIMESTEPPER_NODE_INTERIOR_HPP_

#include "DESolver_TimeTreeNode_Base.hpp"


namespace sweet
{

/*
 * Helper class for interior tree node
 *
 * It provides default member variables which are often used
 */
class DESolver_TimeTreeNode_NodeInteriorHelper	:
		public DESolver_TimeTreeNode_Base
{
protected:
	double _timestep_size;
	double &dt = _timestep_size;

	// DE term to evaluate
	std::vector<std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>> _timeTreeNodes;
	std::vector<DESolver_TimeTreeNode_Base::EvalFun> _evalFuns;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _tmpDataContainer;


public:
	DESolver_TimeTreeNode_NodeInteriorHelper()	:
		_timestep_size(-1)
	{
	}

	virtual
	~DESolver_TimeTreeNode_NodeInteriorHelper()
	{
	}

	bool _helperSetupConfigAndGetTimeStepperEval(
		const sweet::DESolver_Config_Base &i_deTermConfig,
		const std::string &i_timeStepperEvalName,
		DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper,

		const std::string &i_evalName = "integration"
	)
	{
		_evalFuns.resize(_timeTreeNodes.size());
		for (std::size_t i = 0; i < _evalFuns.size(); i++)
		{
			_timeTreeNodes[i]->setupConfigAndGetTimeStepperEval(i_deTermConfig, i_evalName, _evalFuns[i]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}

		// default setup
		DESolver_TimeTreeNode_Base::_helperSetupConfigAndGetTimeStepperEval(
				i_timeStepperEvalName,
				o_timeStepper
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	inline
	void setTimeStepSize(double i_dt)	override
	{
		_timestep_size = i_dt;

		for (auto &i : _timeTreeNodes)
		{
			i->setTimeStepSize(i_dt);
		}
	}


	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		for (auto &i : _timeTreeNodes)
		{
			i->shackRegistration(io_shackDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(*i);
		}

		return true;
	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			delete _tmpDataContainer[i];

		_tmpDataContainer.clear();
	}


	/*
	 * Helper function to call the functions provided by pointers
	 *
	 * This calls the 1st one in the vector
	 */
	inline
	void evalTimeStepper(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		(_timeTreeNodes[0].get()->*_evalFuns[0])(i_U, o_U, i_simulationTime);
	}

	/*
	 * Helper function to call the functions provided by pointers
	 */
	inline
	void evalTimeStepper(
			std::size_t i_id,
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		(_timeTreeNodes[i_id].get()->*_evalFuns[i_id])(i_U, o_U, i_simulationTime);
	}

};

}

#endif
