#include "PDESWESphere_na_vd.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_na_vd::PDESWESphere_na_vd()	:
	shackPDESWESphere(nullptr),
	ops(nullptr),
	dt(-1)
{
	setEvalAvailable("tendencies");
}

PDESWESphere_na_vd::~PDESWESphere_na_vd()
{
}


bool PDESWESphere_na_vd::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_na_vd::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("na_vd");
	return retval;

}

std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> PDESWESphere_na_vd::getNewInstance()
{
	return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_na_vd);
}


bool PDESWESphere_na_vd::setupConfigAndGetTimeStepperEval(
		const sweet::DESolver_Config_Base &i_deTermConfig,
		const std::string &i_timeStepperEvalName,
		DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
)
{
	const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

	ops = myConfig.ops;

	// default setup
	DESolver_TimeTreeNode_Base::_helperSetupConfigAndGetTimeStepperEval(
			i_timeStepperEvalName,
			o_timeStepper
		);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

	return true;
}

void PDESWESphere_na_vd::setTimeStepSize(double i_dt)
{
	dt = i_dt;
}

void PDESWESphere_na_vd::clear()
{
}

/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_na_vd::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_time_stamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(ops != nullptr);
	assert(shackPDESWESphere != nullptr);

	sweet::SphereData_Physical U_u_phys, U_v_phys;
	ops->vrtdiv_to_uv(i_U.vrt, i_U.div, U_u_phys, U_v_phys);

	sweet::SphereData_Physical U_div_phys = i_U.div.toPhys();
	o_U.phi_pert = ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U.phi_pert.toPhys());
	o_U.vrt = ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U.vrt.toPhys());
	o_U.div = ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U.div.toPhys());
}
