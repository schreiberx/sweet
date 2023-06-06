#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPLICITRUNGEKUTTA_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPLICITRUNGEKUTTA_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeInteriorHelper.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_ExplicitRungeKutta	:
	public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_ExplicitRungeKutta>
{
private:
	/*
	 * We use an enum to identify different RK implementations
	 */
	enum ERKMethod
	{
		INVALID = -1,
		ERK1 = 1,	// 1st order forward Euler
		ERK2_MIDPOINT,
		ERK2_HEUN,	// Heun's method
		ERK2_RALSTON,	// Ralston's method
		ERK2_RALSTON_CC,	// Ralston's method, CC version
		ERK3,
		ERK4,		// 4th order classical RK

	};
	ERKMethod _rkMethodID;


	// Order of Runge-Kutta method
	int _order;

	// Particular RK method
	std::string _method;

	// Runge-Kutta stage storages
	int _rkNumStages;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _rkStageDataContainer;

public:
	DESolver_TimeStepper_ExplicitRungeKutta()	:
		_rkMethodID(INVALID),
		_order(-1),
		_method("std"),
		_rkNumStages(-1)
	{
		setEvalAvailable("integration");
	}

	~DESolver_TimeStepper_ExplicitRungeKutta()
	{
		clear();
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("erk");
		retval.push_back("ERK");
		retval.push_back("explicitRungeKutta");
		return retval;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 1)
			return error.set("Exactly one time node term needs to be given"+getNewLineDebugMessage());

		_rkMethodID = INVALID;
		_rkNumStages = -1;

		if (_method == "heun")
		{
			if (_order != -1 && _order != 2)
				return error.set("Order of Heun's method must be 2");

			_rkMethodID = ERK2_HEUN;
			_order = 2;
			_rkNumStages = 2;
		}
		else if (_method == "midpoint")
		{
			if (_order != -1 && _order != 2)
				return error.set("Order of Midpoint method must be 2");

			_rkMethodID = ERK2_MIDPOINT;
			_order = 2;
			_rkNumStages = 2;
		}
		else if (_method == "ralston")
		{
			if (_order != -1 && _order != 2)
				return error.set("Order of Ralston's method must be 2");

			_rkMethodID = ERK2_RALSTON;
			_order = 2;
			_rkNumStages = 2;
		}
		else if (_method == "ralstoncc")
		{
			if (_order != -1 && _order != 2)
				return error.set("Order of Ralston's method must be 2");

			_rkMethodID = ERK2_RALSTON_CC;
			_order = 2;
			_rkNumStages = 2;
		}
		else if (_method == "std" || _method == "default")
		{
			if (_order < 1 || _order > 4)
				return error.set("Order of ERK method needs to be 1, 2, 3 or 4"+getNewLineDebugMessage());

			switch(_order)
			{
				case 1:	_rkMethodID = ERK1; _rkNumStages = 1; break;
				case 2:	_rkMethodID = ERK2_MIDPOINT; _rkNumStages = 2; break;
				case 3:	_rkMethodID = ERK3; _rkNumStages = 3; break;
				case 4:	_rkMethodID = ERK4; _rkNumStages = 4; break;
			}
		}
		else
		{
			return error.set("Unknown method '"+_method+"'");
		}
		assert(_rkNumStages >= 1);


		if (_rkMethodID == INVALID)
			return error.set("Invalid time stepping method");

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
				if (_timeTreeNodes.size() != 0)
					return error.set("a 2nd timestepper was provided!"+a->getNewLineDebugMessage());

				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() > 0)
					return error.set("Only one DETerm is supported"+a->getNewLineDebugMessage());

				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "order" || a->key == "o")
				{
					a->getValue(_order);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}

				if (a->key == "method" || a->key == "m")
				{
					_method = a->value;
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
				"tendencies"
			);


		/*
		 * Setup buffers for RK stage solutions
		 */
		_rkStageDataContainer.resize(_rkNumStages);

		for (std::size_t i = 0; i < _rkStageDataContainer.size(); i++)
			_rkStageDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		/*
		 * Setup temporary buffers
		 */
		if (_rkNumStages > 1)
			_tmpDataContainer.resize(1);

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		return true;

	}

	void clear() override
	{
		for (std::size_t i = 0; i < _rkStageDataContainer.size(); i++)
			delete _rkStageDataContainer[i];
		_rkStageDataContainer.clear();

		DESolver_TimeTreeNode_NodeInteriorHelper::clear();
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ExplicitRungeKutta);
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ExplicitRungeKutta(*this));
	}


	bool _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		switch(_rkMethodID)
		{
		case ERK1:	return _eval_timeIntegration_ERK1(i_U, o_U, i_simulationTime);
		case ERK2_MIDPOINT:	return _eval_timeIntegration_ERK2_Midpoint(i_U, o_U, i_simulationTime);
		case ERK2_HEUN:	return _eval_timeIntegration_ERK2_Heun(i_U, o_U, i_simulationTime);
		case ERK2_RALSTON:	return _eval_timeIntegration_ERK2_Ralston(i_U, o_U, i_simulationTime);
		case ERK2_RALSTON_CC:	return _eval_timeIntegration_ERK2_RalstonCC(i_U, o_U, i_simulationTime);
		case ERK3:	return _eval_timeIntegration_ERK3(i_U, o_U, i_simulationTime);
		case ERK4:	return _eval_timeIntegration_ERK4(i_U, o_U, i_simulationTime);
		default: return error.set("Wrong RK evaluation");
		}
	}

private:
	bool _eval_timeIntegration_ERK1(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(
				0,
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestepSize, *_rkStageDataContainer[0]);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}

private:
	bool _eval_timeIntegration_ERK2_Midpoint(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
		// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
		/*
		 * c     a
		 * 0   |
		 * 1/2 | 1/2
		 * --------------
		 *     | 0   1    b
		 */
		double a2[1] = {0.5};
		double b[2] = {0.0, 1.0};
		double c[1] = {0.5};

		// STAGE 1
		evalTimeStepper(
				0,
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
			);

		// STAGE 2
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		evalTimeStepper(
				0,
				*_tmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulationTime + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[1], *_rkStageDataContainer[1]);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}

private:
	bool _eval_timeIntegration_ERK2_Heun(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		// See https://en.wikipedia.org/wiki/Heun%27s_method
		/*
		 * c     a
		 * 0   |
		 * 1   | 1
		 * --------------
		 *     | 1/2   1/2    b
		 */
		double a2[1] = {1.0};
		double b[2] = {0.5, 0.5};
		double c[1] = {1.0};

		// STAGE 1
		evalTimeStepper(
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
			);

		// STAGE 2
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulationTime + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(dt*b[1], *_rkStageDataContainer[1]);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}


private:
	bool _eval_timeIntegration_ERK2_Ralston(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 * http://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
		 *
		 * See table below Heun's method https://en.wikipedia.org/wiki/Heun%27s_method
		 */
		/*
		 * c     a
		 * 0   |
		 * 2/3 | 2/3
		 * --------------
		 *     | 1/4   3/4  b
		 */
		double a2[1] = {2.0/3.0};
		double b[2] = {1.0/4.0, 3.0/4.0};
		double c[1] = {2.0/3.0};

		// STAGE 1
		evalTimeStepper(
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
		);

		// STAGE 2
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulationTime + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(dt*b[1], *_rkStageDataContainer[1]);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}


private:
	bool _eval_timeIntegration_ERK2_RalstonCC(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 * This is not the original Ralston method, but the one from the book
		 * "Numerical Methods for Engineers" by Chapra and Canale
		 *
		 * This is discussed here:
		 * http://sachinashanbhag.blogspot.com/2016/04/ralstons-method-controversy.html
		 */
		/*
		 * c     a
		 * 0   |
		 * 3/4 | 3/4
		 * --------------
		 *     | 1/3   2/3  b
		 */
		double a2[1] = {3.0/4.0};
		double b[2] = {1.0/3.0, 2.0/3.0};
		double c[1] = {3.0/4.0};

		// STAGE 1
		evalTimeStepper(
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
			);

		// STAGE 2
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulationTime + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(dt*b[1], *_rkStageDataContainer[1]);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}

private:
	bool _eval_timeIntegration_ERK3(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
		// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
		/*
		 * c     a
		 * 0   |
		 * 1/3 | 1/3
		 * 2/3 | 0    2/3
		 * --------------
		 *     | 1/4  0   3/4
		 */
		double a2[1] = {1.0/3.0};
		double a3[2] = {0.0, 2.0/3.0};
		double b[3] = {1.0/4.0, 0.0, 3.0/4.0};
		double c[2] = {1.0/3.0, 2.0/3.0};

		// STAGE 1
		evalTimeStepper(
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
			);

		// STAGE 2
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulationTime + c[0]*dt
			);

		// STAGE 3
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a3[1], *_rkStageDataContainer[1]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[2],
				i_simulationTime + c[1]*dt
			);

		// Closure
		o_U.op_setVectorPlusScalarMulVector(i_U, _timestepSize*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(_timestepSize*b[1], *_rkStageDataContainer[1]);
		o_U.op_addScalarMulVector(_timestepSize*b[2], *_rkStageDataContainer[2]);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}

private:
	bool _eval_timeIntegration_ERK4(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
		// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
		/*
		 * c     a
		 * 0   |
		 * 1/2 | 1/2
		 * 1/2 | 0    1/2
		 * 1   | 0    0    1
		 * --------------
		 *     | 1/6  1/3  1/3  1/6
		 */
		double a2[1] = {0.5};
		double a3[2] = {0.0, 0.5};
		double a4[3] = {0.0, 0.0, 1.0};
		double b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
		double c[3] = {0.5, 0.5, 1.0};

		// STAGE 1
		evalTimeStepper(
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
			);

		// STAGE 2
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulationTime + c[0]*dt
			);

		// STAGE 3
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a3[1], *_rkStageDataContainer[1]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[2],
				i_simulationTime + c[1]*dt
			);


		// STAGE 4
		_tmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a4[2], *_rkStageDataContainer[2]
		);

		evalTimeStepper(
				*_tmpDataContainer[0],
				*_rkStageDataContainer[3],
				i_simulationTime + c[2]*dt
			);

		// Closure
		o_U.op_setVectorPlusScalarMulVector(i_U, _timestepSize*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(_timestepSize*b[1], *_rkStageDataContainer[1]);
		o_U.op_addScalarMulVector(_timestepSize*b[2], *_rkStageDataContainer[2]);
		o_U.op_addScalarMulVector(_timestepSize*b[3], *_rkStageDataContainer[3]);

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "ERK(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << newPrefix << "  method: " << _method << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
