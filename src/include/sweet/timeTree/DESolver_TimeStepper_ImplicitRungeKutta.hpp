#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_IMPLICIT_RUNGEKUTTA_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_IMPLICIT_RUNGEKUTTA_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeInteriorHelper.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_ImplicitRungeKutta	:
		public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_ImplicitRungeKutta>
{
private:
	/*
	 * We use an enum to identify different RK implementations
	 */
	enum IRKMethod
	{
		INVALID = -1,
		IRK1 = 1,	// 1st order forward Euler
		IRK2,
	};
	IRKMethod _rkMethodID;


	// Order of Runge-Kutta method
	int _order;

	// Particular RK method
	std::string _method;

	// Runge-Kutta stage storages
	//int _rkNumStages;

	// Damping factor for 2nd order IRK CN method. 0.5 means no damping
	double _crank_nicolson_damping_factor = 0.5;

	double _dt_explicit = -1;
	double _dt_implicit = -1;

public:
	DESolver_TimeStepper_ImplicitRungeKutta()	:
		_rkMethodID(INVALID),
		_order(-1),
		_method("std")
	{
		setEvalAvailable("integration");
	}

	~DESolver_TimeStepper_ImplicitRungeKutta()
	{
		clear();
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("irk");
		retval.push_back("IRK");
		retval.push_back("implicitRungeKutta");
		return retval;
	}


	bool _setupArgumentInternals()
	{
		_rkMethodID = INVALID;

		if (_order != -1)
		{
			if (_order == 1)
			{
				_method = "std";
			}
			else if (_order == 2)
			{
				_method = "crank_nicolson";
			}
			else
			{
				return error.set("Order of IRK method provided, only 1 or 2 supported"+getNewLineDebugMessage());
			}
		}

		if (_method == "cn" || _method == "crank_nicolson")
		{
			if (_order != -1 && _order != 2)
				return error.set("Order of Crank-Nicolson's method must be 2"+getNewLineDebugMessage());

			_rkMethodID = IRK2;
			_order = 2;
		}
		else if (_method == "std")
		{
			if (_order < 1 || _order > 2)
				return error.set("Order of IRK method needs to be 1 or 2"+getNewLineDebugMessage());

			_rkMethodID = IRK1;
			_order = 1;
		}
		else
		{
			return error.set("Unknown method '"+_method+"'");
		}

		if (_rkMethodID == INVALID)
			return error.set("Invalid time stepping method");

		if (_timeTreeNodes.size() != 1)
			return error.set("DE Term not specified for time stepper"+getNewLineDebugMessage());

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
#if 1
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
				if (timestepperArgFound)
					return error.set("a 2nd timestepper was provided!"+a->getNewLineDebugMessage());

				i_tsAssemblation.assembleTimeStepperByFunction(a->function, _timeTreeNodeEulerBackward);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

				timestepperArgFound = true;
				break;
#endif

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

				if (a->key == "damping" || a->key == "cn_damping")
				{
					a->getValue(_crank_nicolson_damping_factor);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
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
		/*
		 * Manually setup the backward Euler
		 */
		_timeTreeNodes.resize(_order);
		_evalFuns.resize(_order);

		_timeTreeNodes[0]->setupConfigAndGetTimeStepperEval(i_deTermConfig, "eulerBackward", _evalFuns[0]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

		/*
		 * Setup temporary buffers
		 */
		if (_order == 2)
		{
			// Create new instance for explicit evaluation
			_timeTreeNodes[1] = _timeTreeNodes[0]->getInstanceCopy();

			// Crank-Nicolson
			_tmpDataContainer.resize(2);

			// Create new instance
			_timeTreeNodes[1]->setupConfigAndGetTimeStepperEval(i_deTermConfig, "tendencies", _evalFuns[1]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[1]);
		}

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		// Return time stepper for this routine
		DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
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
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ImplicitRungeKutta);
	}


	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ImplicitRungeKutta(*this));
	}

	void setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		// Not required for explicit time stepper, but we do it
		if (_order == 1)
		{
			_timeTreeNodes[0]->setTimeStepSize(i_dt);
			_dt_explicit = i_dt;
		}
		else if (_order == 2)
		{
			// Crank-Nicolson: Use damping factor for explicit time stepper
			_dt_explicit = i_dt*(1.0-_crank_nicolson_damping_factor);
			_dt_implicit = i_dt*_crank_nicolson_damping_factor;

			// Implicit one
			_timeTreeNodes[0]->setTimeStepSize(_dt_implicit);

			// Explicit one
			_timeTreeNodes[1]->setTimeStepSize(_dt_explicit);


		}
		else
		{
			SWEETError("Internal error");
		}

	}

	bool _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		switch(_rkMethodID)
		{
		case IRK1:	return _eval_timeIntegration_IRK1(i_U, o_U, i_simulationTime);
		case IRK2:	return _eval_timeIntegration_IRK2(i_U, o_U, i_simulationTime);
		default: return error.set("Internal error: Wrong IRK Method");
		}
	}

private:
	bool _eval_timeIntegration_IRK1(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		return evalTimeStepper(0, i_U, o_U, i_simulationTime);
	}

private:
	/*
	 * Crank-Nicolson method:
	 *
	 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
	 *
	 * With q the CN damping facor with no damping for q=0.5
	 */
	bool _eval_timeIntegration_IRK2(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		// Forward Euler time step
		evalTimeStepper(1, i_U, *_tmpDataContainer[0], i_simulationTime);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		_tmpDataContainer[1]->op_setVectorPlusScalarMulVector(
				i_U,
				_dt_explicit,
				*_tmpDataContainer[0]
			);

		// Backward Euler step
		evalTimeStepper(0, *_tmpDataContainer[1], o_U, i_simulationTime+dt*_crank_nicolson_damping_factor);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "IRK(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << newPrefix << "  method: " << _method << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
