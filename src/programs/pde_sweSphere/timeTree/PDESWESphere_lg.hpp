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
#include <sweet/timeTree/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_CastHelper.hpp>
#include "PDESWESphere_DataContainer.hpp"
#include "PDESWESphere_DESolver_Config.hpp"


class PDESWESphere_lg	:
	public sweet::DESolver_TimeTreeNode_Base,
	public sweet::DESolver_TimeTreeNode_CastHelper<
			PDESWESphere_DataContainer,
			PDESWESphere_DESolver_Config
		>
{
private:
	sweet::ShackSphereDataOps *shackSphereDataOps;
	ShackPDESWESphere *shackPDESWESphere;
	const sweet::SphereOperators *ops;

	double _dt;

	sweet::ExpFunctions<double> expFunction;

public:
	PDESWESphere_lg();
	~PDESWESphere_lg();

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override;

	virtual
	const std::vector<std::string> getNodeNames();

	std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> getNewInstance() override;

	virtual
	bool setupConfig(
		const sweet::DESolver_Config_Base &i_deTermConfig
	) override;

	void setTimeStepSize(double i_dt) override;

	void clear() override;

	/*
	 * Return the time tendencies of the PDE term
	 */
	void eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_time_stamp
	)	override;

	/*
	 * Compute an exponential integration for a given exp term
	 */
	void eval_exponential(
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_time_stamp
	)	override;
};


#endif
