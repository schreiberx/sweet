#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_


/*
 * Generic includes
 */
#include <sweet/core/ErrorBase.hpp>
#include <sweet/expIntegration/ExpFunction.hpp>
#include "../ShackPDESWESphere.hpp"

/*
 * Time tree node related includes
 */
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeLeafHelper.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_CastHelper.hpp>
#include "PDESWESphere_DataContainer.hpp"
//#include "PDESWESphere_DataContainerComplex.hpp"
#include "PDESWESphere_DESolver_Config.hpp"

#include "../time/PDESWESphereTS_lg_exp_direct.hpp"

class PDESWESphere_lg	:
	public sweet::DESolver_TimeTreeNode_NodeLeafHelper<PDESWESphere_lg>,
	public sweet::DESolver_TimeTreeNode_CastHelper<
			PDESWESphere_DataContainer,
			PDESWESphere_DESolver_Config
		>
{
private:
	ShackPDESWESphere *_shackPDESWESphere;
	sweet::ShackSphereDataOps *_shackSphereDataOps;
	const sweet::SphereOperators *_ops;
	const sweet::SphereOperatorsComplex *_opsComplex;

	sweet::ExpFunction<double> _expFunction;

#if 0
	sweet::DESolver_TimeTreeNode_CastHelper<
				PDESWESphere_DataContainerComplex,
				PDESWESphere_DESolver_Config
			> _castComplex;
#endif

	/*!
	 * Complex-valued time step size for complex-valued backward Euler
	 *
	 * This is used for REXI solvers of the form
	 * 	U_1 = (I-dt*L)^{-1} U_0
	 */
	std::complex<double> _rexiTermAlpha;
	std::complex<double> _rexiTermBeta;
	std::complex<double> _rexiTermGamma;
	bool _rexiTermGammaActive;


public:
	PDESWESphere_lg();
	~PDESWESphere_lg();

	PDESWESphere_lg(const PDESWESphere_lg &i_val);

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override;

	virtual
	const std::vector<std::string> getNodeNames() override;

	virtual
	bool setupConfigAndGetTimeStepperEval(
		const sweet::DESolver_Config_Base &i_deTermConfig,
		EVAL_TYPES i_evalType,
		DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
	) override;

	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override;

	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
	) override;

	void clear() override;


	/*
	 * Return the time tendencies of the PDE term
	 */
public:
	bool _eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	) override;

	/*
	 * Return the backward Euler time step
	 */
public:
	bool _eval_eulerBackward(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	) override;

	/*
	 * Return evaluation of backward Euler time step with complex data
	 */
	bool _eval_rexiTerm(
			const sweet::DESolver_DataContainer_Base &i_U_,
			sweet::DESolver_DataContainer_Base &o_U_,
			double i_timeStamp
	) override;

	/*
	 * Compute an exponential integration for a given exp term
	 */
private:
	bool _eval_exponential(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	) override;
};


#endif
