#include "PDESWESphere_nr_uv.hpp"
#include <vector>
#include <sweet/core/sphere/SphereOperators.hpp>


PDESWESphere_nr_uv::PDESWESphere_nr_uv()	:
	_shackPDESWESphere(nullptr),
	_ops(nullptr)
{
	setEvalAvailable(EVAL_TENDENCIES);
}


PDESWESphere_nr_uv::~PDESWESphere_nr_uv()
{
}

PDESWESphere_nr_uv::PDESWESphere_nr_uv(
		const PDESWESphere_nr_uv &i_value
)	:
	DESolver_TimeTreeNode_NodeLeafHelper(i_value)
{
	_shackPDESWESphere = i_value._shackPDESWESphere;
	_ops = i_value._ops;
}

bool PDESWESphere_nr_uv::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	_shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	return true;
}


const std::vector<std::string> PDESWESphere_nr_uv::getNodeNames()
{
	std::vector<std::string> retval;
	retval.push_back("nr");
	retval.push_back("nr_uv");
	return retval;

}


bool PDESWESphere_nr_uv::setupConfigAndGetTimeStepperEval(
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

void PDESWESphere_nr_uv::clear()
{
	DESolver_TimeTreeNode_NodeLeafHelper::clear();
}

/*
 * Return the time tendencies of the PDE term
 */
bool PDESWESphere_nr_uv::_eval_tendencies(
		const sweet::DESolver_DataContainer_Base &i_U_,
		sweet::DESolver_DataContainer_Base &o_U_,
		double i_timeStamp
)
{
	const PDESWESphere_DataContainer &i_U = cast(i_U_);
	PDESWESphere_DataContainer &o_U = cast(o_U_);

	assert(_ops != nullptr);
	assert(_shackPDESWESphere != nullptr);

	o_U.phi_pert = -sweet::SphereData_Spectral(i_U.phi_pert.toPhys()*i_U.div.toPhys());
	o_U.vrt.spectral_set_zero();
	o_U.div.spectral_set_zero();

	return true;
}
