#include "PDESWESphere_na_uv.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_na_uv::PDESWESphere_na_uv()	:
	shackPDESWESphere(nullptr),
	ops(nullptr),
	dt(-1)
{
	setEvalAvailable("tendencies");
}

PDESWESphere_na_uv::~PDESWESphere_na_uv()
{
}


bool PDESWESphere_na_uv::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_na_uv::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("na");
	retval.push_back("na_uv");
	return retval;

}

std::shared_ptr<sweet::DESolver_TimeTreeNode_Base> PDESWESphere_na_uv::getNewInstance()
{
	return std::shared_ptr<sweet::DESolver_TimeTreeNode_Base>(new PDESWESphere_na_uv);
}


bool PDESWESphere_na_uv::setupConfig(
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

void PDESWESphere_na_uv::setTimeStepSize(double i_dt)
{
	dt = i_dt;
}

void PDESWESphere_na_uv::clear()
{
}

/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_na_uv::eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_time_stamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	sweet::SphereData_Physical U_u_phys, U_v_phys;
	ops->vrtdiv_to_uv(i_U.vrt, i_U.div, U_u_phys, U_v_phys);



	sweet::SphereData_Physical U_div_phys = i_U.div.toPhys();
	o_U.phi_pert = ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, i_U.phi_pert.toPhys());

	/*
	 * Velocity
	 */
	sweet::SphereData_Physical vrtg = i_U.vrt.toPhys();

	sweet::SphereData_Physical u_nl = U_u_phys*vrtg;
	sweet::SphereData_Physical v_nl = U_v_phys*vrtg;

	sweet::SphereData_Spectral vrt, div;
	ops->uv_to_vrtdiv(u_nl, v_nl, vrt, div);
	o_U.vrt = div;

	o_U.div = vrt;
	o_U.div -= ops->laplace(0.5*(U_u_phys*U_u_phys+U_v_phys*U_v_phys));
}
