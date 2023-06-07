#include "PDESWESphere_lc.hpp"
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_lc::PDESWESphere_lc()	:
	_shackPDESWESphere(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}

PDESWESphere_lc::~PDESWESphere_lc()
{
}


PDESWESphere_lc::PDESWESphere_lc(
		const PDESWESphere_lc &i_value
)	:
	DESolver_TimeTreeNode_NodeLeafHelper(i_value)
{
	_shackPDESWESphere = i_value._shackPDESWESphere;
	_ops = i_value._ops;
}


bool PDESWESphere_lc::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	_shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_lc::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("lc");
	return retval;

}


bool PDESWESphere_lc::setupConfigAndGetTimeStepperEval(
	const sweet::DESolver_Config_Base &i_deTermConfig,
	EVAL_TYPES i_evalType,
	DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
)
{
	const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;

	if (_shackPDESWESphere->sphere_use_fsphere)
		fg = _ops->getFG_fSphere(_shackPDESWESphere->sphere_fsphere_f0);
	else
		fg = _ops->getFG_rotatingSphere(_shackPDESWESphere->sphere_rotating_coriolis_omega);

	ug.setup(_ops->sphereDataConfig);
	vg.setup(_ops->sphereDataConfig);

	// default setup
	DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
			i_evalType,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}


void PDESWESphere_lc::clear()
{
	DESolver_TimeTreeNode_NodeLeafHelper::clear();
}


/*
 * Return the time tendencies of the PDE term
 */
bool PDESWESphere_lc::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	/*
	 * step 1a
	 */
	_ops->vrtdiv_to_uv(i_U.vrt, i_U.div, ug, vg);

	/*
	 * step 1b
	 */
	sweet::SphereData_Physical tmp_u = ug*fg;
	sweet::SphereData_Physical tmp_v = vg*fg;

	/*
	 * step 1c
	 */
	_ops->uv_to_vrtdiv(tmp_u, tmp_v, o_U.div, o_U.vrt);

	/*
	 * step 1d
	 */
	o_U.vrt *= -1.0;


	/*
	 * step 1e
	 * Nothing to do
	 */

	/*
	 * step 2a
	 * Zero tendencies
	 */
	o_U.phi_pert.spectral_set_zero();

	return true;
}
