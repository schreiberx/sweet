#include "PDESWESphere_lg.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_lg::PDESWESphere_lg()	:
	shackPDESWESphere(nullptr),
	shackSphereDataOps(nullptr),
	ops(nullptr)
{
	setEvalAvailable("tendencies");
	setEvalAvailable("eulerBackward");
	setEvalAvailable("exponential");
}

PDESWESphere_lg::~PDESWESphere_lg()
{
}


PDESWESphere_lg::PDESWESphere_lg(
		const PDESWESphere_lg &i_val
)	:
	DESolver_TimeTreeNode_NodeLeafHelper(*this)
{
	shackSphereDataOps = i_val.shackSphereDataOps;
	shackPDESWESphere = i_val.shackPDESWESphere;
	ops = i_val.ops;

	setEvalAvailable("tendencies");
	setEvalAvailable("eulerBackward");
	setEvalAvailable("exponential");
}


bool PDESWESphere_lg::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	shackSphereDataOps = io_shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_lg::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("lg");
	return retval;

}


bool PDESWESphere_lg::setupConfigAndGetTimeStepperEval(
	const sweet::DESolver_Config_Base &i_deTermConfig,
	const std::string &i_timeStepperEvalName,
	DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
)
{
	const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

	ops = myConfig.ops;

	// default setup
	DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
			i_timeStepperEvalName,
			o_timeStepper
		);

	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	if (i_timeStepperEvalName == "exponential")
		expFunction.setup("phi0");

	return true;
}


void PDESWESphere_lg::clear()
{
	DESolver_TimeTreeNode_NodeLeafHelper::clear();
}


/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_lg::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	assert(ops != nullptr);
	assert(shackPDESWESphere != nullptr);

	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

	// TODO: Write in a way which directly writes output to output array
	o_U.phi_pert = -gh*i_U.div;
	o_U.div = -ops->laplace(i_U.phi_pert);
	o_U.vrt.spectral_set_zero();
}



/*
 * Evaluate exponential of linear term
 */
void PDESWESphere_lg::_eval_eulerBackward(
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
	double GH = shackPDESWESphere->h0*shackPDESWESphere->gravitation;


	sweet::SphereData_Spectral rhs = i_U.div + ops->implicit_L(i_U.phi_pert, _dt);
	o_U.div = ops->implicit_helmholtz(rhs, -GH*_dt*_dt, shackSphereDataOps->sphere_radius);
	o_U.phi_pert = i_U.phi_pert - _dt*GH*o_U.div;
	o_U.vrt = i_U.vrt;
}



/*
 * Evaluate exponential of linear term
 */
void PDESWESphere_lg::_eval_exponential(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(expFunction.isSetup());

	/*
	 * Using exponential integrators, we must compute an
	 */
	double ir = 1.0/shackSphereDataOps->sphere_radius;

	/*
	 * Average geopotential
	 */
	double GH = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

	/*
	 * See doc/time_integration/exponential_integration/swe_sphere_direct_exp_int_for_L_nonrotating_sphere/
	 */

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= ops->sphereDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = ops->sphereDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= ops->sphereDataConfig->spectral_modes_n_max; n++)
		{
			double D = (double)n*((double)n+1.0)*ir*ir;

			const std::complex<double> &i_phi = i_U.phi_pert.spectral_space_data[idx];
			const std::complex<double> &i_div = i_U.div.spectral_space_data[idx];
			const std::complex<double> &i_vrt = i_U.vrt.spectral_space_data[idx];

			std::complex<double> &o_phi = o_U.phi_pert.spectral_space_data[idx];
			std::complex<double> &o_div = o_U.div.spectral_space_data[idx];
			std::complex<double> &o_vrt = o_U.vrt.spectral_space_data[idx];

			if (D == 0)
			{
				o_phi = i_phi;
				o_div = i_div;
				o_vrt = i_vrt;

				idx++;
				continue;
			}


			// result will be imaginary only!
			std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(-D*GH));

			std::complex<double> l0 = -sqrt_DG/(-2*GH) * i_phi + 0.5*i_div;
			std::complex<double> l1 = +sqrt_DG/(-2*GH) * i_phi + 0.5*i_div;

			std::complex<double> tmp;
			expFunction.eval(-_dt*sqrt_DG, tmp);
			l0 *= tmp;

			expFunction.eval(_dt*sqrt_DG, tmp);
			l1 *= tmp;

			o_phi = -GH/sqrt_DG * (l1 - l0);
			o_div = l0 + l1;
			o_vrt = i_vrt;

			idx++;
		}
	}
}
