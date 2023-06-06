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
template <typename MainInteriorNodeClass>
class DESolver_TimeTreeNode_NodeInteriorHelper	:
		public DESolver_TimeTreeNode_Base
{
protected:
	double _timestepSize;
	double &dt = _timestepSize;

	// DE term to evaluate
	std::vector<std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>> _timeTreeNodes;
	std::vector<DESolver_TimeTreeNode_Base::EvalFun> _evalFuns;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _tmpDataContainer;


public:
	DESolver_TimeTreeNode_NodeInteriorHelper()	:
		_timestepSize(-1)
	{
	}

	virtual
	~DESolver_TimeTreeNode_NodeInteriorHelper()
	{
		clear();
	}


	/*
	 * Copy constructor
	 *
	 * Simply copy the raw data over here
	 */
	DESolver_TimeTreeNode_NodeInteriorHelper(
		const DESolver_TimeTreeNode_NodeInteriorHelper &i_src
	)
	{
		_timestepSize = i_src._timestepSize;

		// registered evaluation functions
		_registeredEvalFunctions = i_src._registeredEvalFunctions;

		// Data containers
		_tmpDataContainer.resize(_tmpDataContainer.size());
		for (std::size_t i = 0; i < i_src._tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_src._tmpDataContainer[i]->getInstanceNew();

		// Time tree nodes
		assert(i_src._timeTreeNodes.size() != 0);
		_timeTreeNodes.resize(i_src._timeTreeNodes.size());
		for (std::size_t i = 0; i < _timeTreeNodes.size(); i++)
			_timeTreeNodes[i] = i_src._timeTreeNodes[i]->getInstanceCopy();

		// _evalFuns is handled later (hopefully)
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
		DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
				i_timeStepperEvalName,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	inline
	void setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		for (auto &i : _timeTreeNodes)
		{
			i->setTimeStepSize(_timestepSize);
		}
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
	bool evalTimeStepper(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		assert(_timeTreeNodes[0] != nullptr);
		assert(_evalFuns[0] != nullptr);

		(_timeTreeNodes[0].get()->*_evalFuns[0])(i_U, o_U, i_simulationTime);
#if SWEET_DEBUG
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);
#endif
		return true;
	}

	/*
	 * Helper function to call the functions provided by pointers
	 */
	inline
	bool evalTimeStepper(
			std::size_t i_id,
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		assert(_timeTreeNodes.size() > i_id);
		assert(_evalFuns.size() > i_id);

		(_timeTreeNodes[i_id].get()->*_evalFuns[i_id])(i_U, o_U, i_simulationTime);
#if SWEET_DEBUG
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i_id]);
#endif
		return true;
	}


	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new MainInteriorNodeClass);
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy()	override
	{
		MainInteriorNodeClass *m = new MainInteriorNodeClass(static_cast<MainInteriorNodeClass&>(*this));
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(m);
	}
};

}

#endif
