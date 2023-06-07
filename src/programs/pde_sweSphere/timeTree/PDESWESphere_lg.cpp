#include "PDESWESphere_lg.hpp"


#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>

PDESWESphere_lg::PDESWESphere_lg()	:
	_shackPDESWESphere(nullptr),
	_shackSphereDataOps(nullptr),
	_ops(nullptr),
	_opsComplex(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
	setEvalAvailable(EVAL_EULER_BACKWARD);
	setEvalAvailable(EVAL_EXPONENTIAL);
	setEvalAvailable(EVAL_REXI_TERM);
}

PDESWESphere_lg::~PDESWESphere_lg()
{
}


PDESWESphere_lg::PDESWESphere_lg(
		const PDESWESphere_lg &i_value
)	:
	DESolver_TimeTreeNode_NodeLeafHelper(i_value),
	_rexiTermAlpha(-999999),
	_rexiTermBeta(-666666)
{
	_shackSphereDataOps = i_value._shackSphereDataOps;
	_shackPDESWESphere = i_value._shackPDESWESphere;
	_ops = i_value._ops;
	_opsComplex = i_value._opsComplex;
}


bool PDESWESphere_lg::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	_shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	_shackSphereDataOps = io_shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}

const std::vector<std::string> PDESWESphere_lg::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("lg");
	return retval;

}


bool PDESWESphere_lg::setupByKeyValue(
		const std::string &i_key,
		const std::string &i_value
)
{
	if (i_key == "expIntegrationFunction")
	{
		if (_expFunction.functionName != "")
			return error.set("Function name for expFunction is already set ('"+_expFunction.functionName+"')");

		if (i_value == "")
			return error.set("Empty string for expFunction given");

		_expFunction.setup(i_value);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(_expFunction);

		return true;
	}

	if (i_key == "rexiTermPreallocation")
	{
		return true;
#if 0
		if (i_value == "true")
		{
			_rexiTermPreallocation = true;
			return true;
		}

		if (i_value == "false")
		{
			_rexiTermPreallocation = false;
			return true;
		}
#endif
		return error.set("Invalid key-value combination of '"+i_key+"' => '"+i_value+"'");
	}

	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


/*!
 * Setup a key which value is complex valued
 */
bool PDESWESphere_lg::setupByKeyValue(
		const std::string &i_key,
		const std::complex<double> &i_value
)
{
	if (i_key == "rexiTermAlpha")
	{
		_rexiTermAlpha = i_value;
		return true;
	}

	if (i_key == "rexiTermBeta")
	{
		_rexiTermBeta = i_value;
		return true;
	}

	return error.set("setupByKeyValue key '"+i_key+"' not supported");
}


bool PDESWESphere_lg::setupConfigAndGetTimeStepperEval(
	const sweet::DESolver_Config_Base &i_deTermConfig,
	EVAL_TYPES i_evalType,
	DESolver_TimeTreeNode_Base::EvalFun &o_timeStepper
)
{
	const PDESWESphere_DESolver_Config& myConfig = cast(i_deTermConfig);

	_ops = myConfig.ops;
	_opsComplex = myConfig.opsComplex;

	assert(_ops != nullptr);
	assert(_opsComplex != nullptr);

	// default setup
	DESolver_TimeTreeNode_Base::_helperGetTimeStepperEval(
			i_evalType,
			o_timeStepper
		);

	ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

	if (i_evalType == EVAL_EXPONENTIAL)
	{
		if (_expFunction.functionName == "")
		{
			// set default to phi0
			_expFunction.setup("phi0");
		}
	}

	return true;
}


void PDESWESphere_lg::clear()
{
	DESolver_TimeTreeNode_NodeLeafHelper::clear();
}


/*
 * Return the time tendencies of the PDE term
 */
bool PDESWESphere_lg::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	assert(_ops != nullptr);
	assert(_shackPDESWESphere != nullptr);

	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	double gh = _shackPDESWESphere->gravitation * _shackPDESWESphere->h0;

	// TODO: Write in a way which directly writes output to output array
	o_U.phi_pert = -gh*i_U.div;
	o_U.div = -_ops->laplace(i_U.phi_pert);
	o_U.vrt.spectral_set_zero();

	return true;
}



/*
 * Evaluate exponential of linear term
 */
bool PDESWESphere_lg::_eval_eulerBackward(
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
	double GH = _shackPDESWESphere->h0*_shackPDESWESphere->gravitation;


	sweet::SphereData_Spectral rhs = i_U.div + _ops->implicit_L(i_U.phi_pert, _dt);
	o_U.div = _ops->implicit_helmholtz(rhs, -GH*_dt*_dt, _shackSphereDataOps->sphere_radius);
	o_U.phi_pert = i_U.phi_pert - _dt*GH*o_U.div;
	o_U.vrt = i_U.vrt;

	return true;
}



/*!
 * Evaluate exponential of linear term
 */
bool PDESWESphere_lg::_eval_exponential(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(_expFunction.isSetup());

	/*
	 * Using exponential integrators, we must compute an
	 */
	double ir = 1.0/_shackSphereDataOps->sphere_radius;

	/*
	 * Average geopotential
	 */
	double GH = _shackPDESWESphere->gravitation*_shackPDESWESphere->h0;

	/*
	 * See doc/time_integration/exponential_integration/swe_sphere_direct_exp_int_for_L_nonrotating_sphere/
	 */

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= _ops->sphereDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = _ops->sphereDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= _ops->sphereDataConfig->spectral_modes_n_max; n++)
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
			_expFunction.eval(-_dt*sqrt_DG, tmp);
			l0 *= tmp;

			_expFunction.eval(_dt*sqrt_DG, tmp);
			l1 *= tmp;

			o_phi = -GH/sqrt_DG * (l1 - l0);
			o_div = l0 + l1;
			o_vrt = i_vrt;

			idx++;
		}
	}

	return true;
}



/*!
 * Evaluate REXI term
 *
 * o_U = \beta (\alpha + dt*L)^-1 i_U
 */
bool PDESWESphere_lg::_eval_rexiTerm(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

#if 0
	const std::complex<double> &alpha = _rexiTermAlpha;
	const std::complex<double> &beta = _rexiTermBeta;

	sweet::SphereData_SpectralComplex i_phi_pert_cplx = sweet::Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_U.phi_pert);
	sweet::SphereData_SpectralComplex i_vrt_cplx = sweet::Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_U.vrt);
	sweet::SphereData_SpectralComplex i_div_cplx = sweet::Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_U.div);

	/*
	 * Rearrange problem to
	 *
	 * o_U = \beta/\alpha (I + dt/\alpha*L)^-1 i_U
	 */
	const std::complex<double> dtComplex = _dt/alpha;

	// avg. geopotential
	double GH = _shackPDESWESphere->h0*_shackPDESWESphere->gravitation;

	sweet::SphereData_SpectralComplex rhs = i_div_cplx + _opsComplex->implicit_L(i_phi_pert_cplx, dtComplex);
	sweet::SphereData_SpectralComplex o_div_cplx = _opsComplex->implicit_helmholtz(rhs, -GH*dtComplex*dtComplex, _shackSphereDataOps->sphere_radius);
	sweet::SphereData_SpectralComplex o_phi_pert_cplx = i_phi_pert_cplx - dtComplex*GH*o_div_cplx;
	sweet::SphereData_SpectralComplex o_vrt_cplx = i_vrt_cplx;

	o_phi_pert_cplx *= beta/alpha;
	o_div_cplx *= beta/alpha;
	o_vrt_cplx *= beta/alpha;

	o_U.phi_pert = sweet::Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(o_phi_pert_cplx);
	o_U.div = sweet::Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(o_div_cplx);
	o_U.vrt = sweet::Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(o_vrt_cplx);

#else

	sweet::SphereData_SpectralComplex phi1(i_U.phi_pert.sphereDataConfig);
	sweet::SphereData_SpectralComplex vrt1(i_U.phi_pert.sphereDataConfig);
	sweet::SphereData_SpectralComplex div1(i_U.phi_pert.sphereDataConfig);

	sweet::SphereData_SpectralComplex phi0 = sweet::Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_U.phi_pert);
	sweet::SphereData_SpectralComplex vrt0 = sweet::Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_U.vrt);
	sweet::SphereData_SpectralComplex div0 = sweet::Convert_SphereDataSpectral_To_SphereDataSpectralComplex::physical_convert(i_U.div);

	/*
	 * Preprocessing:
	 * U0* = U0 * beta/alpha
	 * dt_implicit = -timestep_size/alpha
	 */
	phi0 *= _rexiTermBeta/_rexiTermAlpha;
	div0 *= _rexiTermBeta/_rexiTermAlpha;
	vrt0 *= _rexiTermBeta/_rexiTermAlpha;

	double gh0 = _shackPDESWESphere->h0*_shackPDESWESphere->gravitation;

#if 1
	std::complex<double> dt_implicit = -_dt/_rexiTermAlpha;

	sweet::SphereData_SpectralComplex rhs = div0 + _opsComplex->implicit_L(phi0, dt_implicit);
	div1 = _opsComplex->implicit_helmholtz(rhs, gh0*dt_implicit*dt_implicit, _shackSphereDataOps->sphere_radius);

	phi1 = phi0 - dt_implicit*gh0*div1;
	vrt1 = vrt0;

#else

	std::complex<double> dt_two_omega = dt_implicit*0.0;

	sweet::SphereData_SpectralComplex rhs = div0 + opsComplex.implicit_FJinv(vrt0, dt_two_omega) + opsComplex.implicit_L(phi0, dt_implicit);
	div1 = sphSolverComplexDiv.solve(rhs);
	phi1 = phi0 - dt_implicit*gh0*div1;
	vrt1 = opsComplex.implicit_Jinv(vrt0 - opsComplex.implicit_F(div1, dt_two_omega), dt_two_omega);

#endif

	o_U.phi_pert = sweet::Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(phi1);
	o_U.vrt = sweet::Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(vrt1);
	o_U.div = sweet::Convert_SphereDataSpectralComplex_To_SphereDataSpectral::physical_convert_real(div1);

#endif
	return true;
}


