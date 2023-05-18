#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LC_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LC_HPP_


#include <sweet/timeTree/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/timeTree/DESolver_TimeTreeNode_BaseHelper.hpp>
#include <sweet/timeTree/DESolver_DataContainer_Base.hpp>
#include "PDESWESphere_DataContainer.hpp"
#include "PDESWESphere_DESolver_Config.hpp"

#include "../ShackPDESWESphere.hpp"
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>


class PDESWESphere_lc :
	public sweet::DESolver_TimeTreeNode_Base,
	public sweet::DESolver_TimeTreeNode_BaseHelper<
			PDESWESphere_lc,
			PDESWESphere_DataContainer,
			PDESWESphere_DESolver_Config
		>
{
private:
	ShackPDESWESphere *shackPDESWESphere;
	const sweet::SphereOperators *ops;

	double dt;
	sweet::SphereData_Physical fg;

	sweet::SphereData_Physical ug;
	sweet::SphereData_Physical vg;

public:
	PDESWESphere_lc()	:
		shackPDESWESphere(nullptr),
		ops(nullptr),
		dt(-1)
	{
	}

	~PDESWESphere_lc()	override
	{
	}

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)	override
	{
		shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}

public:
	const std::vector<std::string> getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("lc");
		return retval;

	}

	std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> getNewInstance() override
	{
		return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_lc);
	}


	virtual
	bool setupConfig(
		const sweet::DESolver_Config_Base &i_deTermConfig
	) override
	{
		const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

 		ops = myConfig.ops;

		if (shackPDESWESphere->sphere_use_fsphere)
			fg = ops->getFG_fSphere(shackPDESWESphere->sphere_fsphere_f0);
		else
			fg = ops->getFG_rotatingSphere(shackPDESWESphere->sphere_rotating_coriolis_omega);

		ug.setup(ops->sphereDataConfig);
		vg.setup(ops->sphereDataConfig);

		return true;
	}

	void setTimeStepSize(double i_dt) override
	{
		dt = i_dt;
	}


	void clear() override
	{
	}


	/*
	 * Return the time tendencies of the PDE term
	 */
	void eval_tendencies(
			const sweet::DESolver_DataContainer_Base &i_U,
			sweet::DESolver_DataContainer_Base &o_U,
			double i_time_stamp
	)	override
	{
		/*
		 * step 1a
		 */
		ops->vrtdiv_to_uv(cast(i_U).vrt, cast(i_U).div, ug, vg);

		/*
		 * step 1b
		 */
		sweet::SphereData_Physical tmp_u = ug*fg;
		sweet::SphereData_Physical tmp_v = vg*fg;

		/*
		 * step 1c
		 */
		ops->uv_to_vrtdiv(tmp_u, tmp_v, cast(o_U).div, cast(o_U).vrt);

		/*
		 * step 1d
		 */
		cast(o_U).vrt *= -1.0;


		/*
		 * step 1e
		 * Nothing to do
		 */

		/*
		 * step 2a
		 * Zero tendencies
		 */
		cast(o_U).phi_pert.spectral_set_zero();
	}
};


#endif
