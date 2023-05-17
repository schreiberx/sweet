#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPLICITRUNGEKUTTA_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_EXPLICITRUNGEKUTTA_HPP_

#include <vector>
#include <string>
#include "PDESolver_DataContainer_Base.hpp"
#include "PDESolver_TimeStepper_Base.hpp"
#include "sweet/core/shacks/ShackDictionary.hpp"

class MyTimeStepper_ExplicitRungeKutta	:
		public sweet::PDESolver_TimeStepper_Base
{
private:
	double _timestep_size;
	double &dt = _timestep_size;

	/*
	 * We use an enum to identify different RK implementations
	 */
	enum RKMethod
	{
		INVALID = -1,
		ERK1 = 1,	// 1st order forward Euler
		ERK2,
		ERK3,
		ERK4		// 4th order classical RK
	};
	RKMethod _rkMethodID;

	// PDE term to evaluate
	std::shared_ptr<PDESolver_TimeStepper_Base> _timeStepper;
	std::shared_ptr<sweet::PDESolver_PDETerm_Base> _pdeTerm;

	// Order of Runge-Kutta method
	int _order;

	// Runge-Kutta stage storages
	int _rkNumStages;

	// Number of stages to allocate buffers
	std::vector<sweet::PDESolver_DataContainer_Base*> _rkStageDataContainer;
	std::vector<sweet::PDESolver_DataContainer_Base*> _rkTmpDataContainer;


public:
	MyTimeStepper_ExplicitRungeKutta()	:
		_timestep_size(-1),
		_rkMethodID(INVALID),
		_order(-1),
		_rkNumStages(-1)
	{
	}

	~MyTimeStepper_ExplicitRungeKutta()
	{
		clear();
	}

	const std::vector<std::string>
	getImplementedTimeSteppers()	override
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
		return true;
	}

	bool setupWithArguments(
			std::shared_ptr<sweet::PDESolver_PDETerm_Base> &i_pdeTerm,
			int i_order
	)
	{
		_order = i_order;
		_pdeTerm = i_pdeTerm;

		_rkMethodID = INVALID;
		_rkNumStages = -1;

		switch(_order)
		{
			case 1:	_rkMethodID = ERK1; _rkNumStages = 1; break;
			case 2:	_rkMethodID = ERK2; _rkNumStages = 2; break;
			case 3:	_rkMethodID = ERK3; _rkNumStages = 3; break;
			case 4:	_rkMethodID = ERK4; _rkNumStages = 4; break;
		}

		if (_rkMethodID == INVALID)
			return error.set("Invalid time stepping method");

		return true;
	}

	virtual
	bool setupFunction(
			std::shared_ptr<sweet::PDESolver_TimeStepping_Tree::Function> &i_function,
			sweet::PDESolver_TimeStepping_Assemblation &i_tsAssemblation
	)
	{
		bool timestepperArgFound = false;
		//bool pdeTermFound = false;

		int i = 0;
		for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
		{
			sweet::PDESolver_TimeStepping_Tree::Argument *a = (*iter).get();

			switch(a->argType)
			{

			case sweet::PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_FUNCTION:
				if (a->key != "fun")
					return error.set("Only key 'fun' supported for a function\n"+a->debug_message);
				// continue with ARG_TYPE_FUNCTION

			case sweet::PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_FUNCTION:
				if (timestepperArgFound)
					return error.set("a 2nd timestepper was provided!\n"+a->debug_message);

				i_tsAssemblation.assembleTimeStepperByFunction(a->function, _timeStepper);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

				timestepperArgFound = true;
				break;

			case sweet::PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "order" || a->key == "o")
				{
					a->getValue(_order);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
				}
				else
				{
					return error.set("Key not supported\n"+a->debug_message);
				}
				break;

			case sweet::PDESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				i_tsAssemblation.assemblePDETermByString(
						a->value,
						_pdeTerm
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
				break;

			default:
				return error.set("Internal error");
			}


			i++;
		}

		if (_order <= 0 || _order >= 4)
			return error.set("Order of ERK method needs to be 1, 2, 3 or 4");

		return true;
	}

	bool setupDataContainers(
			const sweet::PDESolver_DataContainer_Base &i_U
	)	override
	{
		clear();

		if (_rkNumStages <= 0 || _rkNumStages > 4)
			return error.set("Invalid order for RK time stepping");

		/*
		 * Setup buffers for RK stage solutions
		 */
		_rkStageDataContainer.resize(_rkNumStages);

		for (std::size_t i = 0; i < _rkStageDataContainer.size(); i++)
			_rkStageDataContainer[i] = i_U.getNewInstance();

		if (_rkNumStages > 1)
			_rkTmpDataContainer.resize(1);

		/*
		 * Setup temporary buffers
		 */
		for (std::size_t i = 0; i < _rkTmpDataContainer.size(); i++)
			_rkTmpDataContainer[i] = i_U.getNewInstance();

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

	std::shared_ptr<PDESolver_TimeStepper_Base> getNewInstance()	override
	{
		return std::shared_ptr<PDESolver_TimeStepper_Base>(new MyTimeStepper_ExplicitRungeKutta);
	}

	void setTimestepSize(double i_dt)	override
	{
		_timestep_size = i_dt;
	}

	void eval_timeIntegration(
			sweet::PDESolver_DataContainer_Base &i_U,
			sweet::PDESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)	override
	{
		switch(_rkMethodID)
		{
		case ERK1:	_eval_timeIntegration_ERK1(i_U, o_U, i_simulation_time);	return;
		case ERK2:	_eval_timeIntegration_ERK2(i_U, o_U, i_simulation_time);	return;
		case ERK3:	_eval_timeIntegration_ERK3(i_U, o_U, i_simulation_time);	return;
		case ERK4:	_eval_timeIntegration_ERK4(i_U, o_U, i_simulation_time);	return;
		default: ;
		}
	}

private:
	void _eval_timeIntegration_ERK1(
			sweet::PDESolver_DataContainer_Base &i_U,
			sweet::PDESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_pdeTerm->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);
	}

private:
	void _eval_timeIntegration_ERK2(
			sweet::PDESolver_DataContainer_Base &i_U,
			sweet::PDESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_pdeTerm->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

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
		_pdeTerm->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		// STAGE 2
		_rkTmpDataContainer[0]->op_setVectorPlusScalarMulVector(
				i_U, dt*a2[0], *_rkStageDataContainer[0]
		);

		_pdeTerm->eval_tendencies(
				*_rkTmpDataContainer[0],
				*_rkStageDataContainer[1],
				i_simulation_time + c[0]*dt
			);

		// We simply swap the content of o_U and i_U to add the stage values
		o_U.swap(i_U);

		o_U.op_addScalarMulVector(dt*b[1], *_rkStageDataContainer[1]);
	}

private:
	void _eval_timeIntegration_ERK3(
			sweet::PDESolver_DataContainer_Base &i_U,
			sweet::PDESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_pdeTerm->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);
	}

private:
	void _eval_timeIntegration_ERK4(
			sweet::PDESolver_DataContainer_Base &i_U,
			sweet::PDESolver_DataContainer_Base &o_U,
			double i_simulation_time
	)
	{
		_pdeTerm->eval_tendencies(i_U, *_rkStageDataContainer[0], i_simulation_time);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestep_size, *_rkStageDataContainer[0]);
	}
};



#endif
