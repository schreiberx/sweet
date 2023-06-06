#ifndef SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_ETDnRK_HPP_
#define SRC_PROGRAMS_SIMDATA_MYTIMESTEPPER_ETDnRK_HPP_

#include <vector>
#include <string>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeInteriorHelper.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


namespace sweet
{

class DESolver_TimeStepper_ETDRK	:
	public DESolver_TimeTreeNode_NodeInteriorHelper<DESolver_TimeStepper_ETDRK>
{
private:
	// Order of Runge-Kutta method
	int _order;

	// Number of different phi variants we require (phi0, phi1, ...)
	int _numPhiVariants;

public:
	DESolver_TimeStepper_ETDRK()	:
		_order(-1),
		_numPhiVariants(-1)
	{
		setEvalAvailable("integration");
	}

	~DESolver_TimeStepper_ETDRK()
	{
		clear();
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("etdrk");
		retval.push_back("ETDRK");
		return retval;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 2)
			return error.set("We require two terms (linear + nonlinear) for this time stepper!"+getNewLineDebugMessage());

		if (_order != 1 && _order == 2 && _order == 4)
			return error.set("Only order 1, 2, and 4 supported for ETDRK");

		if (_order == 1)
		{
			_numPhiVariants = 2;
		}
		else if (_order == 2)
		{
			_numPhiVariants = 3;
		}
		else if (_order == 4)
		{
			_numPhiVariants = 7;
		}

		_timeTreeNodes[0].swap(_timeTreeNodes[1]);

		// We abuse the _timeTreeNodes to store also the different variants
		_timeTreeNodes.resize(1+_numPhiVariants);

		/*
		 * _timeTreeNodes[0]:	 Nonlinear part
		 * _timeTreeNodes[1]:	 phi0(...*L)
		 * _timeTreeNodes[2]:	 phi1(...*L)
		 * _timeTreeNodes[3]:	 phi2(...*L)
		 * _timeTreeNodes[4]:	 phi0(...*L)
		 * _timeTreeNodes[5]:	 ups1(...*L)
		 * _timeTreeNodes[6]:	 ups2(...*L)
		 * _timeTreeNodes[7]:	 ups3(...*L)
		 */

		for (int i = 2; i < 1+_numPhiVariants; i++)
		{
			_timeTreeNodes[i] = _timeTreeNodes[1]->getInstanceCopy();
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}

		if (_order == 1)
		{
			assert(_timeTreeNodes.size() == 3);
			_timeTreeNodes[1]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("expIntegrationFunction", "phi1");
		}
		else if (_order == 2)
		{
			assert(_timeTreeNodes.size() == 4);
			_timeTreeNodes[1]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("expIntegrationFunction", "phi1");
			_timeTreeNodes[3]->setupByKeyValue("expIntegrationFunction", "phi2");
		}
		else if (_order == 4)
		{
			assert(_timeTreeNodes.size() == 8);
			_timeTreeNodes[1]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("expIntegrationFunction", "phi1");
			_timeTreeNodes[3]->setupByKeyValue("expIntegrationFunction", "phi2");

			_timeTreeNodes[4]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[5]->setupByKeyValue("expIntegrationFunction", "ups1");
			_timeTreeNodes[6]->setupByKeyValue("expIntegrationFunction", "ups2");
			_timeTreeNodes[7]->setupByKeyValue("expIntegrationFunction", "ups3");
		}

		for (auto &i : _timeTreeNodes)
		{
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i);
		}

		return true;
	}

	inline
	void setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		_timeTreeNodes[0]->setTimeStepSize(_timestepSize);

		if (_order == 1)
		{
			_timeTreeNodes[1]->setTimeStepSize(dt);
			_timeTreeNodes[2]->setTimeStepSize(dt);
			return;
		}
		else if (_order == 2)
		{
			_timeTreeNodes[1]->setTimeStepSize(dt);
			_timeTreeNodes[2]->setTimeStepSize(dt);
			_timeTreeNodes[3]->setTimeStepSize(dt);
			return;
		}
		else if (_order == 4)
		{
			_timeTreeNodes[1]->setTimeStepSize(0.5*dt);
			_timeTreeNodes[2]->setTimeStepSize(0.5*dt);
			_timeTreeNodes[3]->setTimeStepSize(0.5*dt);

			_timeTreeNodes[4]->setTimeStepSize(dt);
			_timeTreeNodes[5]->setTimeStepSize(dt);
			_timeTreeNodes[6]->setTimeStepSize(dt);
			_timeTreeNodes[7]->setTimeStepSize(dt);
			return;
		}

		SWEETError("Internal error");
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
				if (_timeTreeNodes.size() >= 2)
					return error.set("A 3rd timestepper was provided, but only 2 allowed!"+a->getNewLineDebugMessage());

				_timeTreeNodes.push_back(std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::DESolver_TimeStepping_Tree::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() >= 2)
					return error.set("A 3rd timestepper was provided, but only 2 allowed!"+a->getNewLineDebugMessage());

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
		_evalFuns.resize(_timeTreeNodes.size());

		_timeTreeNodes[0]->setupConfigAndGetTimeStepperEval(i_deTermConfig, "tendencies", _evalFuns[0]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

		for (std::size_t i = 1; i < _timeTreeNodes.size(); i++)
		{
			_timeTreeNodes[i]->setupConfigAndGetTimeStepperEval(i_deTermConfig, "exponential", _evalFuns[i]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}

		// default setup
		DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
				i_timeStepperEvalName,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		/*
		 * Setup temporary buffers
		 */
		if (_order == 1)
			_tmpDataContainer.resize(3);
		else if (_order == 2)
			_tmpDataContainer.resize(6);
		else if (_order == 4)
			_tmpDataContainer.resize(16);

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		return true;

	}

	void clear() override
	{
		DESolver_TimeTreeNode_NodeInteriorHelper::clear();
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceNew()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ETDRK);
	}

	std::shared_ptr<DESolver_TimeTreeNode_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<DESolver_TimeTreeNode_Base>(new DESolver_TimeStepper_ETDRK(*this));
	}


	bool _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		bool retval;
		switch(_order)
		{
		case 1:	retval = _eval_timeIntegration_ETDRK1(i_U, o_U, i_simulationTime);	break;
		case 2:	retval = _eval_timeIntegration_ETDRK2(i_U, o_U, i_simulationTime);	break;
		case 4:	retval = _eval_timeIntegration_ETDRK4(i_U, o_U, i_simulationTime);	break;
		default: SWEETError("Internal error");
		}
		return retval;
	}

private:
	inline
	void _evalNL(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(0, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi0(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(1, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi1(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(2, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi2(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(3, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps0(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(4, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps1(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(5, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps2(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(6, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps3(
			const DESolver_DataContainer_Base &i_U,
			DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(7, i_U, o_U, i_simulationTime);
	}


private:
	bool _eval_timeIntegration_ETDRK1(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 *         NNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 *         ============================                               NNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 *         ============================           NNNNNNNNNNNNNNNNNNN =========
		 */
		sweet::DESolver_DataContainer_Base &phi1_FU = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1_FU, i_simulationTime);

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 * NNNNN   ============================           =================== =========
		 */
		o_U.op_setVectorPlusScalarMulVector(phi0_U, _timestepSize, phi1_FU);

		return true;
	}


private:
	bool _eval_timeIntegration_ETDRK2(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 *         NNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 *         ========================                                NNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 *         ========================            NNNNNNNNNNNNNNNNNNN ========
		 */
		sweet::DESolver_DataContainer_Base &phi1_FU = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1_FU, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 * NNNNN   ========================            =================== ========
		 */
		sweet::DESolver_DataContainer_Base &A = *_tmpDataContainer[3];
		A.op_setVectorPlusScalarMulVector(phi0_U, dt, phi1_FU);

		/*
		 * U_{n+1} = A_{n} + \Delta t \psi_{2}(\Delta tL) [ (F(A_{n},t_{n}+\Delta t)-F(U_{n}) ]
		 *           =====                                   NNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &FA = *_tmpDataContainer[4];
		_evalNL(A, FA, i_simulationTime + dt);

		/*
		 * U_{n+1} = A_{n} + \Delta t \psi_{2}(\Delta tL) [ (F(A_{n},t_{n}+\Delta t)-F(U_{n}) ]
		 *           =====            NNNNNNNNNNNNNNNNNNN =====================================
		 */
		FA.op_subVector(FU);
		sweet::DESolver_DataContainer_Base &phi2_X = *_tmpDataContainer[5];
		_evalPhi2(FA, phi2_X, i_simulationTime);

		/*
		 * U_{n+1} = A_{n} + \Delta t \psi_{2}(\Delta tL) [ (F(A_{n},t_{n}+\Delta t)-F(U_{n}) ]
		 * =======   =====            =================== =====================================
		 */
		o_U.op_setVectorPlusScalarMulVector(A, dt, phi2_X);

		return true;
	}


private:
	bool _eval_timeIntegration_ETDRK4(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		double dt_half = dt*0.5;

		/*
		 * Precompute commonly used terms
		 */

		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 *         NNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 *         ============================                                   NNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 *         ============================           NNNNNNNNNNNNNNNNNNNNNNN ========
		 */
		sweet::DESolver_DataContainer_Base &phi1 = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1, i_simulationTime);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 * NNNNN   ============================               ======================= ========
		 */
		sweet::DESolver_DataContainer_Base &A = *_tmpDataContainer[3];
		A.op_setVectorPlusScalarMulVector(phi0_U, dt_half, phi1);

		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 *         ============================                                       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &FA = *_tmpDataContainer[4];
		_evalNL(A, FA, i_simulationTime+dt_half);

		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 *         ============================               ======================= ==============================
		 */
		_evalPhi1(FA, phi1, i_simulationTime);

		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 * =====   ============================               ======================= ==============================
		 */
		sweet::DESolver_DataContainer_Base &B = *_tmpDataContainer[5];
		B.op_setVectorPlusScalarMulVector(phi0_U, dt_half, phi1);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         NNNNNNNNNNNNNNNNNNNNNNNNNNNN                                                                          ==============
		 */
		sweet::DESolver_DataContainer_Base &phi0_A = *_tmpDataContainer[6];
		_evalPhi0(A, phi0_A, i_simulationTime);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         ============================                                            NNNNNNNNNNNNNNNNNNNNNNNNNNNNN ==============
		 */
		sweet::DESolver_DataContainer_Base &FB = *_tmpDataContainer[7];
		_evalNL(B, FB, i_simulationTime+dt_half);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         ============================                                        NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::DESolver_DataContainer_Base &tmp1 = *_tmpDataContainer[8];
#if 1
		tmp1.op_setVector(FB);
		tmp1.op_mulScalar(2.0);
		tmp1.op_subVector(FU);
#else
		tmp1.op_setVectorPlusScalarMulVector(FU, -2.0, FB);
		tmp1.op_mulScalar(-1);
#endif


		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         ============================               NNNNNNNNNNNNNNNNNNNNNNNN =================================================
		 */
		_evalPhi1(tmp1, phi1, i_simulationTime);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 * NNNNN   ============================               ======================== =================================================
		 */
		sweet::DESolver_DataContainer_Base &C = *_tmpDataContainer[9];
		C.op_setVectorPlusScalarMulVector(phi0_A, dt_half, phi1);


		/*
		 * R0 - R3
		 */
		sweet::DESolver_DataContainer_Base &FC = *_tmpDataContainer[10];
		_evalNL(C, FC, i_simulationTime+dt);

		const sweet::DESolver_DataContainer_Base &R0 = i_U;
		sweet::DESolver_DataContainer_Base &R1 = FU;

		sweet::DESolver_DataContainer_Base &R2 = *_tmpDataContainer[11];
		R2.op_setVectorPlusScalarMulVector(FA, 1.0, FB);

		sweet::DESolver_DataContainer_Base &R3 = FC;

		/*
		 * U_{n+1} =
		 * 		\psi_{0}(\Delta tL)R_{0}
		 * 			+ \Delta t
		 * 			(
		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
		 * 				  \upsilon_{3}(\Delta tL) R_{3}
		 * 			)
		 */

		sweet::DESolver_DataContainer_Base &R0_ = *_tmpDataContainer[12];
		_evalUps0(R0, R0_, i_simulationTime);

		sweet::DESolver_DataContainer_Base &R1_ = *_tmpDataContainer[13];
		_evalUps1(R1, R1_, i_simulationTime);

		sweet::DESolver_DataContainer_Base &R2_ = *_tmpDataContainer[14];
		_evalUps2(R2, R2_, i_simulationTime);

		sweet::DESolver_DataContainer_Base &R3_ = *_tmpDataContainer[15];
		_evalUps3(R3, R3_, i_simulationTime);


		o_U.op_setVectorPlusScalarMulVector(R0_, dt, R1_);
		o_U.op_addScalarMulVector(2.0*dt, R2_);
		o_U.op_addScalarMulVector(dt, R3_);

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "ERK(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}

#endif
