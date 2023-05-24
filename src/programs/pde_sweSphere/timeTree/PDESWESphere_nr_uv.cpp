#include "PDESWESphere_nr_uv.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_nr_uv::PDESWESphere_nr_uv()	:
	shackPDESWESphere(nullptr),
	ops(nullptr),
	dt(-1)
{
}

PDESWESphere_nr_uv::~PDESWESphere_nr_uv()
{
}


bool PDESWESphere_nr_uv::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_nr_uv::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("nr");
	retval.push_back("nr_uv");
	return retval;

}

std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> PDESWESphere_nr_uv::getNewInstance()
{
	return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_nr_uv);
}


bool PDESWESphere_nr_uv::setupConfig(
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

void PDESWESphere_nr_uv::setTimeStepSize(double i_dt)
{
	dt = i_dt;
}

void PDESWESphere_nr_uv::clear()
{
}

/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_nr_uv::eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_time_stamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(ops != nullptr);
	assert(shackPDESWESphere != nullptr);

	o_U.phi_pert = sweet::SphereData_Spectral(i_U.phi_pert.toPhys()*i_U.div.toPhys());
	o_U.vrt.spectral_set_zero();
	o_U.div.spectral_set_zero();
}
