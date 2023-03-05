/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SPHERE_ADVECTION_TIMESTEPPERS_HPP_
#define SRC_SPHERE_ADVECTION_TIMESTEPPERS_HPP_


#include "time/PDEAdvectionSphereTS_BaseInterface.hpp"
#include <sweet/core/ErrorBase.hpp>


/**
 * SWE Plane time steppers
 */
class PDEAdvectionSphereTimeSteppers
{
public:
	sweet::ErrorBase error;
	PDEAdvectionSphereTS_BaseInterface *master = nullptr;

private:
	std::vector<PDEAdvectionSphereTS_BaseInterface*> _registered_integrators;


	void _integratorsRegisterAll(
			sweet::ShackDictionary *i_shackDict,
			sweet::SphereOperators *i_op
		);

	void _integratorsFreeAll(
			PDEAdvectionSphereTS_BaseInterface *skip_this = nullptr
		);

public:
	PDEAdvectionSphereTimeSteppers();

	void printImplementedTimesteppingMethods(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);

	bool setup(
			const std::string &i_timestepping_method,
			sweet::ShackDictionary *i_shackDict,
			sweet::SphereOperators *io_ops
	);

	void clear();


	bool shackRegistration(
		sweet::ShackDictionary *io_shackDict
	);

	~PDEAdvectionSphereTimeSteppers();
};


#endif
