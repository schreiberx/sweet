#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LC_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LC_HPP_


#include "PDESolver_PDETerm_Base.hpp"
#include "PDESolver_DataContainer_Base.hpp"
#include "MyDataContainer.hpp"

#include "../pde_sweSphere/ShackPDESWESphere.hpp"
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>


class MyPDETerm_lc	:
		public sweet::PDESolver_PDETerm_Base
{
public:
	sweet::ErrorBase error;

private:
	ShackPDESWESphere *shackPDESWESphere;
	sweet::SphereOperators *ops;

	double dt;
	sweet::SphereData_Physical fg;

public:
	MyPDETerm_lc()	:
		shackPDESWESphere(nullptr),
		ops(nullptr),
		dt(-1)
	{
	}

	~MyPDETerm_lc()	override
	{
	}


private:
	static inline
	MyDataContainer& cast(sweet::PDESolver_DataContainer_Base &i_U)
	{
		return static_cast<MyDataContainer&>(i_U);
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

	const char* getImplementedPDETerm()	override
	{
		return "lc";
	}

	std::shared_ptr<sweet::PDESolver_PDETerm_Base> getNewInstance() override
	{
		return std::shared_ptr<sweet::PDESolver_PDETerm_Base>(new MyPDETerm_lg);
	}

	void setup(
		sweet::SphereOperators *io_ops,
		const sweet::PDESolver_DataContainer_Base &i_u
	) override
	{
		ops = io_ops;

		if (shackPDESWESphere->sphere_use_fsphere)
			fg = ops->getFG_fSphere(shackPDESWESphere->sphere_fsphere_f0);
		else
			fg = ops->getFG_rotatingSphere(shackPDESWESphere->sphere_rotating_coriolis_omega);
	}

	void setTimestepSize(double i_dt) override
	{
		dt = i_dt;
	}

	/*
	 * Return the time tendencies of the PDE term
	 */
	void eval_tendencies(
			sweet::PDESolver_DataContainer_Base &i_U,
			sweet::PDESolver_DataContainer_Base &o_U,
			double i_time_stamp
	)	override
	{
		// TODO: Move to setup
		sweet::SphereData_Physical ug(cast(i_U).phi.sphereDataConfig);
		sweet::SphereData_Physical vg(cast(i_U).phi.sphereDataConfig);

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
		cast(o_U).phi.spectral_set_zero();
	}

	/*
	 * Return the forward Euler evaluation of the term:
	 *
	 * (U^{n+1}-U^{n})/Dt = dU^{n}/dt
	 * <=> U^{n+1} = U^n + Dt* (dU^{n}n/dt)
	 */
#if 0
	void eval_euler_forward(
			PDESolver_DataContainer_Base &i_u,
			PDESolver_DataContainer_Base &o_u
		)	override
	{
	}
#endif
};


#endif /* SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_ */
