/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDEAdvectionSphereTS_na_erk.hpp"




bool PDEAdvectionSphereTS_na_erk::testImplementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	return i_timestepping_method == "na_erk";
}

std::string PDEAdvectionSphereTS_na_erk::getStringId()
{
	return "na_erk";
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDEAdvectionSphereTS_na_erk::euler_timestep_update(
		const sweet::SphereData_Spectral &i_prognostic_field,	///< prognostic variables
		sweet::SphereData_Physical &io_u,
		sweet::SphereData_Physical &io_v,

		sweet::SphereData_Spectral &o_prognostic_field,	///< time updates

		double i_simulation_timestamp
)
{

	/**
	 * We simply compute
	 * 	-DIV(rho*U) = -rho DIV(U) - U.GRAD(rho) = - U.GRAD(rho)
	 * which is the Lagrangian contribution only.
	 *
	 * This is the case because the velocity field is divergence free!!!
	 */
	sweet::SphereData_Spectral phi = i_prognostic_field;

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (shackPDEAdvBenchmark->getVelocities)
	{
		shackPDEAdvBenchmark->getVelocities(io_u, io_v, i_simulation_timestamp, shackPDEAdvBenchmark->getVelocitiesUserData);

		sweet::SphereData_Spectral vrt(phi.sphereDataConfig);
		sweet::SphereData_Spectral div(phi.sphereDataConfig);
		ops->uv_to_vrtdiv(io_u, io_v, vrt, div);
	}

	sweet::SphereData_Physical phig = phi.toPhys();

	sweet::SphereData_Physical tmpg1 = io_u*phig;
	sweet::SphereData_Physical tmpg2 = io_v*phig;

	o_prognostic_field = -ops->uv_to_div(tmpg1, tmpg2);
}



void PDEAdvectionSphereTS_na_erk::run_timestep(
		std::vector<sweet::SphereData_Spectral> &io_prognostic_fields,	///< prognostic variables
		sweet::SphereData_Physical &io_u,
		sweet::SphereData_Physical &io_v,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	for (std::size_t i = 0; i < io_prognostic_fields.size(); i++)
	{
		// standard time stepping
		timestepping_rk.run_timestep_na(
				this,
				&PDEAdvectionSphereTS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
				io_prognostic_fields[i], io_u, io_v,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
}


bool PDEAdvectionSphereTS_na_erk::setup(
	sweet::SphereOperators *io_ops
)
{
	PDEAdvectionSphereTS_BaseInterface::setup(io_ops);

	assert(shackPDEAdvectionTimeDisc != nullptr);
	timestepping_order = shackPDEAdvectionTimeDisc->timestepping_order;

	if (timestepping_order < 1 || timestepping_order > 4)
		return error.set("Invalid time stepping order");

	return true;
}


void PDEAdvectionSphereTS_na_erk::printImplementedTimesteppingMethods(
		std::ostream &o_ostream,
		const std::string &i_prefix
)
{
	o_ostream << i_prefix << " + PDEAdvectionSphereTS_na_erk:" << std::endl;
	o_ostream << i_prefix << "    * 'na_erk'" << std::endl;
}


PDEAdvectionSphereTS_na_erk::~PDEAdvectionSphereTS_na_erk()
{
}

