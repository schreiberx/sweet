#include "PDESWESphere_n.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_n::PDESWESphere_n()	:
	_shackPDESWESphere(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}

PDESWESphere_n::~PDESWESphere_n()
{
}

PDESWESphere_n::PDESWESphere_n(
		const PDESWESphere_n &i_value
)	:
	DESolver_TimeTreeNode_NodeLeafHelper(i_value)
{
	_shackPDESWESphere = i_value._shackPDESWESphere;
	_ops = i_value._ops;
}


bool PDESWESphere_n::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	_shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_n::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("n");
	return retval;

}


bool PDESWESphere_n::setupConfigAndGetTimeStepperEval(
	const sweet::DESolver_Config_Base &i_deTermConfig,
	EVAL_TYPES i_evalType,
	DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
)
{
	const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;

	// default setup
	DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
			i_evalType,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}

void PDESWESphere_n::clear()
{
	DESolver_TimeTreeNode_NodeLeafHelper::clear();
}

/*
 * Return the time tendencies of the PDE term
 */
bool PDESWESphere_n::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(_ops != nullptr);
	assert(_shackPDESWESphere != nullptr);


	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	sweet::SphereData_Physical ug(i_U.phi_pert.sphereDataConfig);
	sweet::SphereData_Physical vg(i_U.phi_pert.sphereDataConfig);

	sweet::SphereData_Physical vrtg = i_U.vrt.toPhys();
	sweet::SphereData_Physical divg = i_U.div.toPhys();
	_ops->vrtdiv_to_uv(i_U.vrt, i_U.div, ug, vg);

	sweet::SphereData_Physical phig = i_U.phi_pert.toPhys();

	sweet::SphereData_Physical tmpg1 = ug*(vrtg/*+fg*/);
	sweet::SphereData_Physical tmpg2 = vg*(vrtg/*+fg*/);

	_ops->uv_to_vrtdiv(tmpg1, tmpg2, o_U.div, o_U.vrt);

	o_U.vrt *= -1.0;

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	sweet::SphereData_Spectral tmpspec(i_U.phi_pert.sphereDataConfig);
	_ops->uv_to_vrtdiv(tmpg1,tmpg2, tmpspec, o_U.phi_pert);

	o_U.phi_pert *= -1.0;

	sweet::SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	tmpspec = /*phig+*/tmpg;

	o_U.div += -_ops->laplace(tmpspec);

	return true;
}
