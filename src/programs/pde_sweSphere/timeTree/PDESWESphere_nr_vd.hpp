#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_NR_VD_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_NR_VD_HPP_


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


class PDESWESphere_nr_vd	:
	public sweet::DESolver_TimeTreeNode_NodeLeafHelper<PDESWESphere_nr_vd>,
	public sweet::DESolver_TimeTreeNode_CastHelper<
			PDESWESphere_DataContainer,
			PDESWESphere_DESolver_Config
		>
{
private:
	ShackPDESWESphere *shackPDESWESphere;
	const sweet::SphereOperators *ops;

public:
	PDESWESphere_nr_vd();
	~PDESWESphere_nr_vd();

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

private:
	void euler_timestep_update_na(
			const sweet::SphereData_Spectral &i_U_phi,
			const sweet::SphereData_Spectral &i_U_vrt,
			const sweet::SphereData_Spectral &i_U_div,

			sweet::SphereData_Spectral &o_phi_t,
			sweet::SphereData_Spectral &o_vrt_t,
			sweet::SphereData_Spectral &o_div_t,

			double i_simulation_timestamp
	);

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
