#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LN_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LN_HPP_


/*
 * Generic includes
 */
#include <sweet/core/ErrorBase.hpp>
#include "../ShackPDESWESphere.hpp"

/*
 * Time tree node related includes
 */
#include <sweet/timeTree/DESolver_TimeTreeNode_NodeLeafHelper.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_CastHelper.hpp>
#include "PDESWESphere_DataContainer.hpp"
#include "PDESWESphere_DESolver_Config.hpp"


class PDESWESphere_ln	:
	public sweet::DESolver_TimeTreeNode_NodeLeafHelper<PDESWESphere_ln>,
	public sweet::DESolver_TimeTreeNode_CastHelper<
			PDESWESphere_DataContainer,
			PDESWESphere_DESolver_Config
		>
{
private:
	ShackPDESWESphere *shackPDESWESphere;
	const sweet::SphereOperators *ops;

	/*
	 * Coriolis effect
	 */
	sweet::SphereData_Physical fg;

	/*
	 * Temporary variables
	 */
	sweet::SphereData_Physical ug;
	sweet::SphereData_Physical vg;

public:
	PDESWESphere_ln();
	~PDESWESphere_ln();

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override;

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
	void _eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	)	override;
};


#endif
