#include "PDESWESphere_ln.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_ln::PDESWESphere_ln()	:
	shackPDESWESphere(nullptr),
	ops(nullptr)
{
	setEvalAvailable("tendencies");
}

PDESWESphere_ln::~PDESWESphere_ln()
{
}


bool PDESWESphere_ln::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_ln::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("ln");
	return retval;

}


bool PDESWESphere_ln::setupConfigAndGetTimeStepperEval(
	const sweet::DESolver_Config_Base &i_deTermConfig,
	const std::string &i_timeStepperEvalName,
	DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
)
{
	const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

	ops = myConfig.ops;

	if (shackPDESWESphere->sphere_use_fsphere)
		fg = ops->getFG_fSphere(shackPDESWESphere->sphere_fsphere_f0);
	else
		fg = ops->getFG_rotatingSphere(shackPDESWESphere->sphere_rotating_coriolis_omega);

	// default setup
	DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
			i_timeStepperEvalName,
			o_timeStepper
		);
	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	return true;
}


void PDESWESphere_ln::clear()
{
	DESolver_TimeTreeNode_NodeLeafHelper::clear();
}

/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_ln::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(ops != nullptr);
	assert(shackPDESWESphere != nullptr);


	double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

	/*
	 * NON-LINEAR SWE
	 *
	 * See
	 * 	Williamson, David L., Drake, John B., Hack, James J., Jakob, Rudiger, & Swarztrauber, Paul N. (1992).
	 * 	A standard test set for numerical approximations to the shallow water equations in spherical geometry.
	 * 	Journal of Computational Physics, 102(1), 211–224. https://doi.org/10.1016/S0021-9991(05)80016-6
	 *
	 * "2.3 Vorticity/Divergence Form"
	 */

	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 */
	sweet::SphereData_Physical phi_pert_phys = i_U.phi_pert.toPhys();

	/*
	 * Step 1a
	 */
	sweet::SphereData_Physical ug, vg;
	ops->vrtdiv_to_uv(i_U.vrt, i_U.div, ug, vg);

	/*
	 * Step 1b
	 */
	sweet::SphereData_Physical vrtg = i_U.vrt.toPhys();

	/*
	 * Step 1c
	 */

	using namespace sweet;

	// left part of eq. (19)
	sweet::SphereData_Physical u_nl = ug*(vrtg+fg);

	// left part of eq. (20)
	sweet::SphereData_Physical v_nl = vg*(vrtg+fg);

	/*
	 * Step 1d
	 */
	// Eq. (21) & left part of Eq. (22)
	ops->uv_to_vrtdiv(u_nl, v_nl, o_U.div, o_U.vrt);


	/*
	 * Step 1e
	 */
	o_U.vrt *= -1.0;

	/*
	 * Step 1f
	 */
	// Right part of Eq. (22)
	sweet::SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	sweet::SphereData_Spectral e = phi_pert_phys+tmpg;

	/*
	 * Step 1g
	 */
	o_U.div -= ops->laplace(e);

	/*
	 * Compute Phi geopotential tendencies
	 */

	/*
	 * Step 2a
	 */
	u_nl = ug*(phi_pert_phys + gh0);
	v_nl = vg*(phi_pert_phys + gh0);

	ops->uv_to_vrtdiv(u_nl,v_nl, e, o_U.phi_pert);

	o_U.phi_pert *= -1.0;
}
