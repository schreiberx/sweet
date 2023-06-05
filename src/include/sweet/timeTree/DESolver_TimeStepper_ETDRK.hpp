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
			_timeTreeNodes[i] = _timeTreeNodes[1]->getInstanceCopy();

		if (_order == 1)
		{
			_timeTreeNodes[1]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("expIntegrationFunction", "phi1");
		}
		else if (_order == 2)
		{
			_timeTreeNodes[1]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("expIntegrationFunction", "phi1");
			_timeTreeNodes[3]->setupByKeyValue("expIntegrationFunction", "phi2");
		}
		else if (_order == 4)
		{
			_timeTreeNodes[1]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("expIntegrationFunction", "phi1");
			_timeTreeNodes[3]->setupByKeyValue("expIntegrationFunction", "phi2");

			_timeTreeNodes[4]->setupByKeyValue("expIntegrationFunction", "phi0");
			_timeTreeNodes[5]->setupByKeyValue("expIntegrationFunction", "ups1");
			_timeTreeNodes[6]->setupByKeyValue("expIntegrationFunction", "ups2");
			_timeTreeNodes[7]->setupByKeyValue("expIntegrationFunction", "ups3");
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
					return error.set("A 3rd timestepper was provided, but only 2 required!"+a->getNewLineDebugMessage());

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


	void _eval_integration(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)	override
	{
		switch(_order)
		{
		case 1:	_eval_timeIntegration_ETDRK1(i_U, o_U, i_simulationTime);	return;
		case 2:	_eval_timeIntegration_ETDRK2(i_U, o_U, i_simulationTime);	return;
		case 4:	_eval_timeIntegration_ETDRK4(i_U, o_U, i_simulationTime);	return;
		default: ;
		}
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
	void _eval_timeIntegration_ETDRK1(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */
#if 0
		sweet::SphereData_Spectral phi0_Un_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_v(sphereDataConfig);
		ts_phi0_rexi.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_fixed_dt,
				i_simulation_timestamp
			);
#endif
		sweet::DESolver_DataContainer_Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

#if 0
		sweet::SphereData_Spectral FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral FUn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);
#endif

		sweet::DESolver_DataContainer_Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

#if 0
		sweet::SphereData_Spectral phi1_FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_v(sphereDataConfig);

		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_fixed_dt,
				i_simulation_timestamp
			);
#endif
		sweet::DESolver_DataContainer_Base &phi1_FU = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1_FU, i_simulationTime);


#if 0
		io_phi_pert = phi0_Un_h + i_fixed_dt*phi1_FUn_h;
		io_vrt = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		io_div = phi0_Un_v + i_fixed_dt*phi1_FUn_v;

		evalTimeStepper(
				0,
				i_U,
				*_rkStageDataContainer[0],
				i_simulationTime
			);

		o_U.op_setVectorPlusScalarMulVector(i_U, _timestepSize, *_rkStageDataContainer[0]);
#endif

		o_U.op_setVectorPlusScalarMulVector(phi0_U, _timestepSize, phi1_FU);
	}


private:
	void _eval_timeIntegration_ETDRK2(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
#if 0
		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		sweet::SphereData_Spectral phi0_Un_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

#endif

		sweet::DESolver_DataContainer_Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

#if 0
		sweet::SphereData_Spectral FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral FUn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

#endif

		sweet::DESolver_DataContainer_Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

#if 0
		sweet::SphereData_Spectral phi1_FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_v(sphereDataConfig);

		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_fixed_dt,
				i_simulation_timestamp
			);
#endif

		sweet::DESolver_DataContainer_Base &phi1_FU = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1_FU, i_simulationTime);

#if 0
		sweet::SphereData_Spectral A_h = phi0_Un_h + i_fixed_dt*phi1_FUn_h;
		sweet::SphereData_Spectral A_u = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		sweet::SphereData_Spectral A_v = phi0_Un_v + i_fixed_dt*phi1_FUn_v;

#endif

		sweet::DESolver_DataContainer_Base &A = *_tmpDataContainer[3];
		A.op_setVectorPlusScalarMulVector(phi0_U, dt, phi1_FU);

#if 0
		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		sweet::SphereData_Spectral FAn_h(sphereDataConfig);
		sweet::SphereData_Spectral FAn_u(sphereDataConfig);
		sweet::SphereData_Spectral FAn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


#endif

		sweet::DESolver_DataContainer_Base &FA = *_tmpDataContainer[4];
		_evalNL(A, FA, i_simulationTime);

#if 0
		sweet::SphereData_Spectral phi2_X_h(sphereDataConfig);
		sweet::SphereData_Spectral phi2_X_u(sphereDataConfig);
		sweet::SphereData_Spectral phi2_X_v(sphereDataConfig);

		ts_phi2_rexi.runTimestep(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h,
				phi2_X_u,
				phi2_X_v,

				i_fixed_dt,
				i_simulation_timestamp
			);

#endif

		FA.op_subVector(FU);
		sweet::DESolver_DataContainer_Base &phi2_X = *_tmpDataContainer[5];
		_evalPhi2(FA, phi2_X, i_simulationTime);

#if 0
		io_phi_pert = A_h + i_fixed_dt*phi2_X_h;
		io_vrt = A_u + i_fixed_dt*phi2_X_u;
		io_div = A_v + i_fixed_dt*phi2_X_v;
#endif
		o_U.op_setVectorPlusScalarMulVector(A, dt, phi2_X);
	}


private:
	void _eval_timeIntegration_ETDRK4(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_simulationTime
	)
	{
		double dt_half = dt*0.5;

		/*
		 * Precompute commonly used terms
		 */
#if 0
		sweet::SphereData_Spectral phi0_Un_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				dt_half,
				i_simulation_timestamp
			);

#endif

		sweet::DESolver_DataContainer_Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

#if 0
		sweet::SphereData_Spectral FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral FUn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);


#endif

		sweet::DESolver_DataContainer_Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * Some commonly shared buffers
		 */
#if 0
		sweet::SphereData_Spectral phi1_h(sphereDataConfig);
		sweet::SphereData_Spectral phi1_u(sphereDataConfig);
		sweet::SphereData_Spectral phi1_v(sphereDataConfig);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 */
		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);
#endif

		sweet::DESolver_DataContainer_Base &phi1 = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1, i_simulationTime);

#if 0
		sweet::SphereData_Spectral A_h = phi0_Un_h + dt_half*phi1_h;
		sweet::SphereData_Spectral A_u = phi0_Un_u + dt_half*phi1_u;
		sweet::SphereData_Spectral A_v = phi0_Un_v + dt_half*phi1_v;
#endif

		sweet::DESolver_DataContainer_Base &A = *_tmpDataContainer[3];
		A.op_setVectorPlusScalarMulVector(phi0_U, dt,  phi1);

		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 */
#if 0
		sweet::SphereData_Spectral FAn_h(sphereDataConfig);
		sweet::SphereData_Spectral FAn_u(sphereDataConfig);
		sweet::SphereData_Spectral FAn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp + dt_half
		);

#endif

		sweet::DESolver_DataContainer_Base &FA = *_tmpDataContainer[4];
		_evalNL(A, FA, i_simulationTime+dt_half);

#if 0
		ts_phi1_rexi.runTimestep(
				FAn_h, FAn_u, FAn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

#endif

		_evalPhi1(FA, phi1, i_simulationTime);

#if 0
		sweet::SphereData_Spectral B_h = phi0_Un_h + dt_half*phi1_h;
		sweet::SphereData_Spectral B_u = phi0_Un_u + dt_half*phi1_u;
		sweet::SphereData_Spectral B_v = phi0_Un_v + dt_half*phi1_v;
#endif

		sweet::DESolver_DataContainer_Base &B = *_tmpDataContainer[5];
		B.op_setVectorPlusScalarMulVector(phi0_U, dt_half, phi1);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 */
#if 0
		sweet::SphereData_Spectral phi0_An_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_An_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_An_v(sphereDataConfig);

		ts_phi0_rexi.runTimestep(
				A_h, A_u, A_v,
				phi0_An_h, phi0_An_u, phi0_An_v,
				dt_half,
				i_simulation_timestamp
			);
#endif

		sweet::DESolver_DataContainer_Base &phi0_A = *_tmpDataContainer[6];
		_evalPhi0(A, phi0_A, i_simulationTime);

#if 0
		sweet::SphereData_Spectral FBn_h(sphereDataConfig);
		sweet::SphereData_Spectral FBn_u(sphereDataConfig);
		sweet::SphereData_Spectral FBn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				B_h, B_u, B_v,
				FBn_h, FBn_u, FBn_v,
				i_simulation_timestamp + dt_half
		);
#endif

		sweet::DESolver_DataContainer_Base &FB = *_tmpDataContainer[7];
		_evalNL(B, FB, i_simulationTime+dt_half);

#if 0
		ts_phi1_rexi.runTimestep(
				2.0*FBn_h - FUn_h,
				2.0*FBn_u - FUn_u,
				2.0*FBn_v - FUn_v,
				phi1_h,	phi1_u,	phi1_v,
				dt_half,
				i_simulation_timestamp
			);
#endif

		sweet::DESolver_DataContainer_Base &tmp1 = *_tmpDataContainer[8];
		tmp1.op_setVectorPlusScalarMulVector(FU, -2.0, FB);
		tmp1.op_mulScalar(-1);

#if 0
		sweet::SphereData_Spectral C_h = phi0_An_h + dt_half*phi1_h;
		sweet::SphereData_Spectral C_u = phi0_An_u + dt_half*phi1_u;
		sweet::SphereData_Spectral C_v = phi0_An_v + dt_half*phi1_v;
#endif

		sweet::DESolver_DataContainer_Base &C = *_tmpDataContainer[9];
		C.op_setVectorPlusScalarMulVector(phi0_A, dt_half,  phi1);

#if 0
		/*
		 * R0 - R3
		 */
		sweet::SphereData_Spectral FCn_h(sphereDataConfig);
		sweet::SphereData_Spectral FCn_u(sphereDataConfig);
		sweet::SphereData_Spectral FCn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				C_h, C_u, C_v,
				FCn_h, FCn_u, FCn_v,
				i_simulation_timestamp + dt
		);
#endif

		sweet::DESolver_DataContainer_Base &FC = *_tmpDataContainer[10];
		_evalNL(C, FC, i_simulationTime+dt);

#if 0

		sweet::SphereData_Spectral R0_h = io_phi_pert;
		sweet::SphereData_Spectral R0_u = io_vrt;
		sweet::SphereData_Spectral R0_v = io_div;

#endif

		const sweet::DESolver_DataContainer_Base &R0 = i_U;

#if 0
		sweet::SphereData_Spectral &R1_h = FUn_h;
		sweet::SphereData_Spectral &R1_u = FUn_u;
		sweet::SphereData_Spectral &R1_v = FUn_v;
#endif

		sweet::DESolver_DataContainer_Base &R1 = FU;

#if 0
		sweet::SphereData_Spectral R2_h = FAn_h + FBn_h;
		sweet::SphereData_Spectral R2_u = FAn_u + FBn_u;
		sweet::SphereData_Spectral R2_v = FAn_v + FBn_v;
#endif

		sweet::DESolver_DataContainer_Base &R2 = *_tmpDataContainer[11];
		R2.op_setVectorPlusScalarMulVector(FA, 1.0, FB);

#if 0
		sweet::SphereData_Spectral &R3_h = FCn_h;
		sweet::SphereData_Spectral &R3_u = FCn_u;
		sweet::SphereData_Spectral &R3_v = FCn_v;
#endif

		sweet::DESolver_DataContainer_Base &R3 = FC;

#if 0
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
		ts_ups0_rexi.runTimestep(
				R0_h, R0_u, R0_v,
				dt,		i_simulation_timestamp
			);

#endif

		sweet::DESolver_DataContainer_Base &R0_ = *_tmpDataContainer[12];
		_evalUps0(R0, R0_, i_simulationTime);

#if 0
		ts_ups1_rexi.runTimestep(
				R1_h, R1_u, R1_v,
				dt,		i_simulation_timestamp
			);
#endif

		sweet::DESolver_DataContainer_Base &R1_ = *_tmpDataContainer[13];
		_evalUps1(R1, R1_, i_simulationTime);

#if 0
		ts_ups2_rexi.runTimestep(
				R2_h, R2_u, R2_v,
				dt,		i_simulation_timestamp
			);
#endif

		sweet::DESolver_DataContainer_Base &R2_ = *_tmpDataContainer[14];
		_evalUps2(R2, R2_, i_simulationTime);

#if 0
		ts_ups3_rexi.runTimestep(
				R3_h, R3_u, R3_v,
				dt,		i_simulation_timestamp
			);
#endif

		sweet::DESolver_DataContainer_Base &R3_ = *_tmpDataContainer[15];
		_evalUps3(R3, R3_, i_simulationTime);

#if 0
		io_phi_pert = R0_h + dt*(R1_h + 2.0*R2_h + R3_h);
		io_vrt = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
		io_div = R0_v + dt*(R1_v + 2.0*R2_v + R3_v);
#endif
		o_U.op_setVectorPlusScalarMulVector(R0_, dt, R1_);
		o_U.op_addScalarMulVector(2.0*dt, R2_);
		o_U.op_addScalarMulVector(dt, R3_);
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
