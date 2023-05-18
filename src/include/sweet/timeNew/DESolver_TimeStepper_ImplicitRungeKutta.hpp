#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_IMPLICIT_RUNGEKUTTA_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_IMPLICIT_RUNGEKUTTA_HPP_

#include <vector>
#include <string>
#include <sweet/timeNew/DESolver_DataContainer_Base.hpp>
#include <sweet/timeNew/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_ImplicitRungeKutta	:
		public sweet::DESolver_TimeTreeNode_Base
{
private:
	double _timestep_size;
	double &dt = _timestep_size;

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

	// DE term to evaluate
	std::shared_ptr<DESolver_TimeTreeNode_Base> _timeTreeNode;

	// Order of Runge-Kutta method
	int _order;

	// Particular RK method
	std::string _method;

	// Runge-Kutta stage storages
	int _rkNumStages;

	// Number of stages to allocate buffers
	std::vector<sweet::DESolver_DataContainer_Base*> _rkTmpDataContainer;

	// Damping factor for 2nd order IRK CN method. 0.5 means no damping
	double _crank_nicolson_damping_factor = 0.5;

public:
	DESolver_TimeStepper_ImplicitRungeKutta()	:
		_timestep_size(-1),
		_rkMethodID(INVALID),
		_order(-1),
		_method("std"),
		_rkNumStages(-1)
	{
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
			return error.set("DE Term not specified for time stepper"+getNewLineDebugMessage());

		_rkMethodID = INVALID;
		_rkNumStages = -1;

		if (_method == "cn" || _method == "crank_nicolson")
		{
			if (_order != -1 && _order != 2)
				return error.set("Order of Crank-Nicolson's method must be 2");

			_rkMethodID = IRK2;
			_order = 2;
			_rkNumStages = 2;
		}
		else if (_method == "std")
		{
			if (_order < 1 || _order > 2)
				return error.set("Order of IRK method needs to be 1 or 2"+getNewLineDebugMessage());

			switch(_order)
			{
				case 1:	_rkMethodID = IRK1; _rkNumStages = 1; break;
				case 2:	_rkMethodID = IRK2; _rkNumStages = 1; break;
			}
		}
		else
		{
			return error.set("Unknown method '"+_method+"'");
		}

		if (!_timeTreeNode->isEvalAvailable("eval_eulerBackward"))
			return error.set("eval_eulerBackward not available");

		if (_order == 2)
			if (!_timeTreeNode->isEvalAvailable("eval_tendencies"))
				return error.set("eval_tendencies not available");


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

				i_tsAssemblation.assembleTimeStepperByFunction(a->function, _timeTreeNode);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

				timestepperArgFound = true;
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

				if (a->key == "damping" || a->key == "cn_damping")
				{
					a->getValue(_crank_nicolson_damping_factor);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNode != nullptr)
					return error.set("Only one DETerm is suppored"+a->getNewLineDebugMessage());

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
		 * Setup temporary buffers
		 */
		if (_order > 1)
			_rkTmpDataContainer.resize(2);

		for (std::size_t i = 0; i < _rkTmpDataContainer.size(); i++)
			_rkTmpDataContainer[i] = i_deConfig.getNewDataContainerInstance();

		return true;

	}

	void clear() override
	{
		for (std::size_t i = 0; i < _rkTmpDataContainer.size(); i++)
			delete _rkTmpDataContainer[i];
		_rkTmpDataContainer.clear();
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getNewInstance()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ImplicitRungeKutta);
	}

	void setTimeStepSize(double i_dt)	override
	{
		_timestep_size = i_dt;

		// Not required for explicit time stepper, but we do it
		if (_order == 1)
			_timeTreeNode->setTimeStepSize(i_dt);
		else if (_order == 2)
		{
			// Crank-Nicolson: Use damping factor for explicit time stepper
			_timeTreeNode->setTimeStepSize(_crank_nicolson_damping_factor*i_dt);
		}
		else
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
		case IRK1:	_eval_timeIntegration_IRK1(i_U, o_U, i_simulation_time);	return;
		case IRK2:	_eval_timeIntegration_IRK2(i_U, o_U, i_simulation_time);	return;
		default: ;
		}
	}

private:
	void _eval_timeIntegration_IRK1(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_timeTreeNode->eval_eulerBackward(i_U, o_U, i_simulation_time);
	}

private:
	/*
	 * Crank-Nicolson method:
	 *
	 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
	 *
	 * With q the CN damping facor with no damping for q=0.5
	 */
	void _eval_timeIntegration_IRK2(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		// Forward Euler step
		_timeTreeNode->eval_tendencies(i_U, *_rkTmpDataContainer[0], i_simulation_time);
		_rkTmpDataContainer[1]->op_setVectorPlusScalarMulVector(
				i_U,
				_crank_nicolson_damping_factor*_timestep_size,
				*_rkTmpDataContainer[0]
			);

		// Backward Euler step
		_timeTreeNode->eval_eulerBackward(*_rkTmpDataContainer[1], o_U, i_simulation_time+dt*_crank_nicolson_damping_factor);
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
