#include "PDESWESphere_l.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_l::PDESWESphere_l()	:
	shackPDESWESphere(nullptr),
	ops(nullptr),
	dt(-1)
{
}

PDESWESphere_l::~PDESWESphere_l()
{
}


bool PDESWESphere_l::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_l::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("l");
	return retval;

}

std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> PDESWESphere_l::getNewInstance()
{
	return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_l);
}


bool PDESWESphere_l::setupConfig(
	const sweet::DESolver_Config_Base &i_deTermConfig
)
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

void PDESWESphere_l::setTimeStepSize(double i_dt)
{
	dt = i_dt;
}

void PDESWESphere_l::clear()
{
}

/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_l::eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_time_stamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(ops != nullptr);
	assert(shackPDESWESphere != nullptr);


	if (!shackPDESWESphere->sphere_use_fsphere)
	{

		double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		/*
		 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
		 */
		/*
		 * Step 1a
		 */
		sweet::SphereData_Physical ug(i_U.phi_pert.sphereDataConfig);
		sweet::SphereData_Physical vg(i_U.phi_pert.sphereDataConfig);
		ops->vrtdiv_to_uv(i_U.vrt, i_U.div, ug, vg);

		/*
		 * Step 1b
		 */
		sweet::SphereData_Physical tmpg1 = ug*fg;
		sweet::SphereData_Physical tmpg2 = vg*fg;

		/*
		 * Step 1c
		 */
		ops->uv_to_vrtdiv(tmpg1, tmpg2, o_U.div, o_U.vrt);

		/*
		 * Step 1d
		 */
		o_U.vrt *= -1.0;

		/*
		 * Step 1e
		 */
		o_U.div += -ops->laplace(i_U.phi_pert);

		/*
		 * DIV on velocity field
		 */
		o_U.phi_pert = (-gh0)*i_U.div;
	}
	else
	{
		double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		o_U.div = -ops->laplace(i_U.phi_pert);

		o_U.vrt = -shackPDESWESphere->sphere_fsphere_f0*i_U.div;
		o_U.div += shackPDESWESphere->sphere_fsphere_f0*i_U.vrt;

		o_U.phi_pert = -gh*i_U.div;
	}
}
