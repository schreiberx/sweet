#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LC_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LC_HPP_


#include <sweet/timeNew/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/timeNew/DESolver_DataContainer_Base.hpp>
#include "PDESWESphere_DataContainer.hpp"

#include "../pde_sweSphere/ShackPDESWESphere.hpp"
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>

#include "PDESWESphere_DESolver_Config.hpp"

class PDESWESphere_lc	:
		public sweet::DESolver_TimeTreeNode_Base
{
private:
	ShackPDESWESphere *shackPDESWESphere;
	const sweet::SphereOperators *ops;

	double dt;
	sweet::SphereData_Physical fg;

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


private:
	static inline
	PDESWESphere_DataContainer& cast(sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<PDESWESphere_DataContainer&>(i_U);
	}

	static inline
	const PDESWESphere_DataContainer& cast(const sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<const PDESWESphere_DataContainer&>(i_U);
	}

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}


	virtual
	const std::vector<std::string> getNodeNames()
	{
		std::vector<std::string> retval;
		retval.push_back("lc");
		return retval;

	}

	std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> getInstanceNew() override
	{
		return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_lc);
	}

	const PDESWESphere_DESolver_Config& cast(
			const sweet::DESolver_Config_Base& i_deTermConfig
	)
	{
		return static_cast<const PDESWESphere_DESolver_Config&>(i_deTermConfig);
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
			double i_timeStamp
	)	override
	{
		// TODO: Move to setup
		sweet::SphereData_Physical ug(cast(i_U).phi_pert.sphereDataConfig);
		sweet::SphereData_Physical vg(cast(i_U).phi_pert.sphereDataConfig);

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

	/*
	 * Return the forward Euler evaluation of the term:
	 *
	 * (U^{n+1}-U^{n})/Dt = dU^{n}/dt
	 * <=> U^{n+1} = U^n + Dt* (dU^{n}n/dt)
	 */
#if 0
	void eval_euler_forward(
			const DESolver_DataContainer_Base &i_u,
			DESolver_DataContainer_Base &o_u
		)	override
	{
	}
#endif
};


#endif /* SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_ */
