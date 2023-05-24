#include "PDESWESphere_lg.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_lg::PDESWESphere_lg()	:
	shackPDESWESphere(nullptr),
	ops(nullptr),
	dt(-1)
{
}

PDESWESphere_lg::~PDESWESphere_lg()
{
}


bool PDESWESphere_lg::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_lg::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("lg");
	return retval;

}

std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> PDESWESphere_lg::getNewInstance()
{
	return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_lg);
}


bool PDESWESphere_lg::setupConfig(
	const sweet::DESolver_Config_Base &i_deTermConfig
)
{
	const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

	ops = myConfig.ops;

	return true;
}

void PDESWESphere_lg::setTimeStepSize(double i_dt)
{
	dt = i_dt;
}

void PDESWESphere_lg::clear()
{
}

/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_lg::eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_u,
		sweet::DESolver_DataContainer_Base &o_u,
		double i_time_stamp
)
{
	assert(ops != nullptr);
	assert(shackPDESWESphere != nullptr);

	const PDESWESphere_DataContainer &i = cast(i_u);
	PDESWESphere_DataContainer &o = cast(o_u);

	double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

	// TODO: Write in a way which directly writes output to output array
	o.phi_pert = -gh*i.div;
	o.div = -ops->laplace(i.phi_pert);
	o.vrt.spectral_set_zero();
}
