#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_L_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_L_HPP_


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
#include "../timeHelpers/SWESphBandedMatrixPhysicalReal.hpp"


class PDESWESphere_l	:
	public sweet::DESolver_TimeTreeNode_NodeLeafHelper<PDESWESphere_l>,
	public sweet::DESolver_TimeTreeNode_CastHelper<
			PDESWESphere_DataContainer,
			PDESWESphere_DESolver_Config
		>
{
private:
	ShackPDESWESphere *shackPDESWESphere;
	sweet::ShackSphereDataOps *shackSphereDataOps;
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

	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolverRealDiv;

public:
	PDESWESphere_l();
	~PDESWESphere_l();

	PDESWESphere_l(
			const PDESWESphere_l &i_val
	);

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

	void clear() override;

	void setTimeStepSize(double i_dt)	override;

	/*
	 * Return the time tendencies of the PDE term
	 */
private:
	bool _eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	)	override;

	/*
	 * Return the time tendencies of the PDE term
	 */
private:
	bool _eval_eulerBackward(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_timeStamp
	)	override;
};


#endif
