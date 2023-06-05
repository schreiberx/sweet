#include "PDESWESphere_l.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_l::PDESWESphere_l()	:
	shackPDESWESphere(nullptr),
	shackSphereDataOps(nullptr),
	ops(nullptr)
{
	setEvalAvailable("tendencies");
	setEvalAvailable("eulerBackward");
}

PDESWESphere_l::~PDESWESphere_l()
{
}


PDESWESphere_l::PDESWESphere_l(
		const PDESWESphere_l &i_val
)	:
	DESolver_TimeTreeNode_NodeLeafHelper(*this)
{
	shackPDESWESphere = i_val.shackPDESWESphere;
	shackSphereDataOps = i_val.shackSphereDataOps;
	ops = i_val.ops;

	setEvalAvailable("tendencies");
	setEvalAvailable("eulerBackward");
}


bool PDESWESphere_l::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	shackSphereDataOps = io_shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_l::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("l");
	return retval;

}


bool PDESWESphere_l::setupConfigAndGetTimeStepperEval(
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

	ug.setup(ops->sphereDataConfig);
	vg.setup(ops->sphereDataConfig);

	if (i_timeStepperEvalName == "eulerBackward")
		if (shackPDESWESphere->sphere_use_fsphere)
			SWEETError("Not supported");


	// default setup
	DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
			i_timeStepperEvalName,
			o_timeStepper
		);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

	return true;
}


void PDESWESphere_l::clear()
{
	DESolver_TimeTreeNode_NodeLeafHelper::clear();
}



void PDESWESphere_l::setTimeStepSize(double i_dt)
{
	_timestepSize = i_dt;

	double two_coriolis = 2.0*shackPDESWESphere->sphere_rotating_coriolis_omega;
	double dt_two_omega = _timestepSize*two_coriolis;
	double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

	sphSolverRealDiv.setup(ops->sphereDataConfig, 4);
	sphSolverRealDiv.solver_component_implicit_J(dt_two_omega);
	sphSolverRealDiv.solver_component_implicit_FJinvF(dt_two_omega);
	sphSolverRealDiv.solver_component_implicit_L(gh0*_timestepSize, _timestepSize, shackSphereDataOps->sphere_radius);
}



/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_l::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
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


/*
 * Evaluate backward Euler timestep
 */
void PDESWESphere_l::_eval_eulerBackward(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	/*
	 * avg. geopotential
	 */
	double gh0 = shackPDESWESphere->h0*shackPDESWESphere->gravitation;

	double dt_two_omega = _dt*2.0*shackPDESWESphere->sphere_rotating_coriolis_omega;

	sweet::SphereData_Spectral rhs = i_U.div + ops->implicit_FJinv(i_U.vrt, dt_two_omega) + ops->implicit_L(i_U.phi_pert, _dt);
	o_U.div = sphSolverRealDiv.solve(rhs);

	o_U.phi_pert = i_U.phi_pert - _dt*gh0*o_U.div;
	o_U.vrt = ops->implicit_Jinv(i_U.vrt - ops->implicit_F(o_U.div, dt_two_omega), dt_two_omega);
}


