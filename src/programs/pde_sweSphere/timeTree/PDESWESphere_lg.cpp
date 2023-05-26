#include "PDESWESphere_lg.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_lg::PDESWESphere_lg()	:
	shackPDESWESphere(nullptr),
	ops(nullptr),
	_dt(-1)
{
	setEvalAvailable("tendencies");
	setEvalAvailable("exponential");
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
	_dt = i_dt;

	expFunction.setup("phi0");
}

void PDESWESphere_lg::clear()
{
}

/*
 * Return the time tendencies of the PDE term
 */
void PDESWESphere_lg::eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_time_stamp
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
void PDESWESphere_lg::eval_exponential(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_time_stamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	/*
	 * Using exponential integrators, we must compute an
	 */
	double ir = 1.0/shackSphereDataOps->sphere_radius;

	/*
	 * avg. geopotential
	 */
	double G = -shackPDESWESphere->h0*shackPDESWESphere->gravitation;

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

			if (D == 0)
			{
				o_U.phi_pert.spectral_space_data[idx] = 0;
				o_U.div.spectral_space_data[idx] = 0;

				idx++;
				continue;
			}

			const std::complex<double> &i_phi = i_U.phi_pert.spectral_space_data[idx];
			const std::complex<double> &i_div = i_U.div.spectral_space_data[idx];

			std::complex<double> &o_phi = o_U.phi_pert.spectral_space_data[idx];
			std::complex<double> &o_div = o_U.div.spectral_space_data[idx];


			// result will be imaginary only!
			std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));

			std::complex<double> l0 = -sqrt_DG/(2*G) * i_phi + 0.5*i_div;
			std::complex<double> l1 = +sqrt_DG/(2*G) * i_phi + 0.5*i_div;

			l0 = expFunction.eval(_dt*(-sqrt_DG))*l0;
			l1 = expFunction.eval(_dt*sqrt_DG)*l1;


			o_phi = -G/sqrt_DG * l0 + G/sqrt_DG* l1;
			o_div = l0 + l1;

			idx++;
		}
	}

	o_U.vrt = i_U.vrt;
}
