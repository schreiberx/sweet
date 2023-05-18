#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPLICITRUNGEKUTTA_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPLICITRUNGEKUTTA_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_ExplicitRungeKutta	:
		public sweet::DESolver_TimeTreeNode_Base
{
private:
	double _timestep_size;
	double &dt = _timestep_size;

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

	// DE term to evaluate
	std::shared_ptr<DESolver_TimeTreeNode_Base> _timeTreeNode;


	// Order of Runge-Kutta method
	int _order;

	// Particular RK method
	std::string _method;

	// Runge-Kutta stage storages
	int _rkNumStages;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _rkStageDataContainer;
	std::vector<sweet::DESolver_DataContainer_Base*> _rkTmpDataContainer;

public:
	DESolver_TimeStepper_ExplicitRungeKutta()	:
		_timestep_size(-1),
		_rkMethodID(INVALID),
		_order(-1),
		_method("std"),
		_rkNumStages(-1)
	{
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

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		_timeTreeNode->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNode);

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNode == nullptr)
			return error.set("Some time node term needs to be given"+getNewLineDebugMessage());

		if (!_timeTreeNode->isEvalAvailable("eval_tendencies"))
			return error.set("eval_tendencies not available in DE term");

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
		else if (_method == "std")
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


		if (_rkMethodID == INVALID)
			return error.set("Invalid time stepping method");

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
			sweet::DESolver_TimeStepping_Tree::Argument *a = (*iter).get();

			switch(a->argType)
			{
#if 0
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION:
				error.set("Time steppers inside this time stepper are not allowed, yet"+a->getNewLineDebugMessage());
				return false;
				break;
#else
			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
				if (a->key != "fun")
					return error.set("Only key 'fun' supported for a function!"+a->getNewLineDebugMessage());
				// continue with ARG_TYPE_FUNCTION

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION:
				if (_timeTreeNode != nullptr)
					return error.set("a 2nd timestepper was provided!"+a->getNewLineDebugMessage());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(a->function, _timeTreeNode);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;
#endif

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "order" || a->key == "o")
				{
					a->getValue(_order);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

				if (a->key == "method" || a->key == "m")
				{
					_method = a->value;
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:

				if (_timeTreeNode != nullptr)
					return error.set("Only one DETerm is supported"+a->getNewLineDebugMessage());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNode
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
		_timeTreeNode->setupConfig(i_deConfig);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNode);

		clear();

		/*
		 * Setup buffers for RK stage solutions
		 */
		_rkStageDataContainer.resize(_rkNumStages);

		for (std::size_t i = 0; i < _rkStageDataContainer.size(); i++)
			_rkStageDataContainer[i] = i_deConfig.getNewDataContainerInstance();

		/*
		 * Setup temporary buffers
		 */
		if (_rkNumStages > 1)
			_rkTmpDataContainer.resize(1);

		for (std::size_t i = 0; i < _rkTmpDataContainer.size(); i++)
			_rkTmpDataContainer[i] = i_deConfig.getNewDataContainerInstance();

		return true;

	}

	void clear() override
	{
		for (std::size_t i = 0; i < _rkStageDataContainer.size(); i++)
			delete _rkStageDataContainer[i];
		_rkStageDataContainer.clear();

		for (std::size_t i = 0; i < _rkTmpDataContainer.size(); i++)
			delete _rkTmpDataContainer[i];
		_rkTmpDataContainer.clear();
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getNewInstance()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ExplicitRungeKutta);
	}

	void setTimeStepSize(double i_dt)	override
	{
		_timestep_size = i_dt;

		// Not required for explicit time stepper, but we do it
		_timeTreeNode->setTimeStepSize(i_dt);
	}

	void eval_timeIntegration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)	override
	{
		switch(_rkMethodID)
		{
		case ERK1:	_eval_timeIntegration_ERK1(i_U, o_U, i_simulation_time);	return;
		case ERK2_MIDPOINT:	_eval_timeIntegration_ERK2_Midpoint(i_U, o_U, i_simulation_time);	return;
		case ERK2_HEUN:	_eval_timeIntegration_ERK2_Heun(i_U, o_U, i_simulation_time);	return;
		case ERK2_RALSTON:	_eval_timeIntegration_ERK2_Ralston(i_U, o_U, i_simulation_time);	return;
		case ERK2_RALSTON_CC:	_eval_timeIntegration_ERK2_RalstonCC(i_U, o_U, i_simulation_time);	return;
		case ERK3:	_eval_timeIntegration_ERK3(i_U, o_U, i_simulation_time);	return;
		case ERK4:	_eval_timeIntegration_ERK4(i_U, o_U, i_simulation_time);	return;
		default: ;
		}
	}

private:
	void _eval_timeIntegration_ERK1(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);
	}

private:
	void _eval_timeIntegration_ERK2_Midpoint(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);


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
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		// STAGE 2
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulation_time + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[1], *_rkStageDataContainer[1]);
	}

private:
	void _eval_timeIntegration_ERK2_Heun(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);


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
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		// STAGE 2
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulation_time + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(dt*b[1], *_rkStageDataContainer[1]);
	}


private:
	void _eval_timeIntegration_ERK2_Ralston(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);


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
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		// STAGE 2
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulation_time + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(dt*b[1], *_rkStageDataContainer[1]);
	}


private:
	void _eval_timeIntegration_ERK2_RalstonCC(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);


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
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		// STAGE 2
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulation_time + c[0]*dt
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, dt*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(dt*b[1], *_rkStageDataContainer[1]);
	}

private:
	void _eval_timeIntegration_ERK3(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
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
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		// STAGE 2
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulation_time + c[0]*dt
			);

		// STAGE 3
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a3[1], *_rkStageDataContainer[1]
		);

		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[2],
				i_simulation_time + c[1]*dt
			);

		// Closure
		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(_timestep_size*b[1], *_rkStageDataContainer[1]);
		o_U.op_addScalarMulVector(_timestep_size*b[2], *_rkStageDataContainer[2]);
	}

private:
	void _eval_timeIntegration_ERK4(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
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
		_timeTreeNode->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		// STAGE 2
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);
		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulation_time + c[0]*dt
			);

		// STAGE 3
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a3[1], *_rkStageDataContainer[1]
		);
		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[2],
				i_simulation_time + c[1]*dt
			);


		// STAGE 4
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a4[2], *_rkStageDataContainer[2]
		);
		_timeTreeNode->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[3],
				i_simulation_time + c[2]*dt
			);

		// Closure
		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size*b[0], *_rkStageDataContainer[0]);
		o_U.op_addScalarMulVector(_timestep_size*b[1], *_rkStageDataContainer[1]);
		o_U.op_addScalarMulVector(_timestep_size*b[2], *_rkStageDataContainer[2]);
		o_U.op_addScalarMulVector(_timestep_size*b[3], *_rkStageDataContainer[3]);
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
