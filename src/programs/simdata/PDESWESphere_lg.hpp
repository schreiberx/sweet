#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_


#include <sweet/timeNew/DESolver_TimeTreeNode_Base.hpp>
#include <sweet/timeNew/DESolver_DataContainer_Base.hpp>
#include <sweet/timeNew/DESolver_Config_Base.hpp>
#include "MyDataContainer.hpp"

#include "../pde_sweSphere/ShackPDESWESphere.hpp"
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>


class PDESWESphere_lg	:
		public sweet::DESolver_TimeTreeNode_Base
{
public:
	sweet::ErrorBase error;

private:
	ShackPDESWESphere *shackPDESWESphere;
	const sweet::SphereOperators *ops;

	double dt;

public:
	PDESWESphere_lg()	:
		shackPDESWESphere(nullptr),
		ops(nullptr),
		dt(-1)
	{
	}

	~PDESWESphere_lg()	override
	{
	}


private:
	static inline
	MyDataContainer& cast(sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<MyDataContainer&>(i_U);
	}

	static inline
	const MyDataContainer& cast(const sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<const MyDataContainer&>(i_U);
	}

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
		retval.push_back("lg");
		return retval;

	}

	std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> getNewInstance() override
	{
		return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_lg);
	}

	virtual
	bool setupConfig(
		const sweet::DESolver_Config_Base &i_deTermConfig
	) override
	{
		//const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

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
			const sweet::DESolver_DataContainer_Base &i_u,
			sweet::DESolver_DataContainer_Base &o_u,
			double i_time_stamp
	)	override
	{
		assert(ops != nullptr);

		const MyDataContainer &i = cast(i_u);
		MyDataContainer &o = cast(o_u);

		double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		// TODO: Write in a way which directly writes output to output array
		o.phi = -gh*i.div;
		o.div = -ops->laplace(i.phi);
		o.vrt.spectral_set_zero();
	}

	/*
	 * Return the forward Euler evaluation of the term:
	 *
	 * (U^{n+1}-U^{n})/Dt = dU^{n}/dt
	 * <=> U^{n+1} = U^n + Dt* (dU^{n}n/dt)
	 */
#if 0
	void eval_euler_forward(
			DESolver_DataContainer_Base &i_u,
			DESolver_DataContainer_Base &o_u
	) override
	{
	}
#endif
};


#endif /* SRC_PROGRAMS_SIMDATA_TIMESTEPPERPDETERM_LG_HPP_ */
