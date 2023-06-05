#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_


/*
 * Generic includes
 */
#include <sweet/expIntegration/ExpFunctions.hpp>
#include <sweet/core/ErrorBase.hpp>
#include "../ShackPDESWESphere.hpp"

/*
 * Time tree node related includes
 */
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeLeafHelper.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_CastHelper.hpp>
#include "PDESWESphere_DataContainer.hpp"
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
	ShackPDESWESphere *shackPDESWESphere;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	const sweet::SphereOperators *ops;

	sweet::ExpFunctions<double> expFunction;

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
		const std::string &i_timeStepperEvalName,
		DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
	) override;

	void clear() override;


	/*
	 * Return the time tendencies of the PDE term
	 */
public:
	void _eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	)	override;

	/*
	 * Return the backward Euler time step
	 */
public:
	void _eval_eulerBackward(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	)	override;

	/*
	 * Compute an exponential integration for a given exp term
	 */
private:
	void _eval_exponential(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	)	override;
};


#endif
