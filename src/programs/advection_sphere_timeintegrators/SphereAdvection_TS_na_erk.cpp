/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../advection_sphere_timeintegrators/SphereAdvection_TS_na_erk.hpp"




bool SphereAdvection_TS_na_erk::implements_timestepping_method(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "na_erk";
}

std::string SphereAdvection_TS_na_erk::string_id()
{
	return "na_erk";
}


void SphereAdvection_TS_na_erk::setup_auto()
{
	setup(simVars.disc.timestepping_order);
}


/*
 * Main routine for method to be used in case of finite differences
 */
void SphereAdvection_TS_na_erk::euler_timestep_update(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vrt,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

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

	SphereData_Spectral phi = i_phi;
	SphereData_Spectral vrt = i_vrt;
	SphereData_Spectral div = i_div;

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (sphereBenchmarks)
	{
		SphereData_Spectral tmp(phi.sphereDataConfig);
		sphereBenchmarks->master->get_reference_state(tmp, vrt, div, i_simulation_timestamp);
	}

	SphereData_Physical ug(phi.sphereDataConfig);
	SphereData_Physical vg(phi.sphereDataConfig);
	op.vortdiv_to_uv(vrt, div, ug, vg);

	SphereData_Physical vrtg = vrt.toPhys();
	SphereData_Physical divg = div.toPhys();

	SphereData_Physical phig = phi.toPhys();

	SphereData_Physical tmpg1 = ug*phig;
	SphereData_Physical tmpg2 = vg*phig;

	SphereData_Spectral tmpspec(phi.sphereDataConfig);
	op.uv_to_vortdiv(tmpg1, tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	o_vort_t.spectral_set_zero();
	o_div_t.spectral_set_zero();
}



void SphereAdvection_TS_na_erk::run_timestep(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const BenchmarksSphereAdvection *i_sphereBenchmarks
)
{
	sphereBenchmarks = i_sphereBenchmarks;

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SphereAdvection_TS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);

	if (simVars.benchmark.getExternalForcesCallback != nullptr)
	{
		// this is just called for cosmetic reasons to update the velocity field
		simVars.benchmark.getExternalForcesCallback(1, simVars.timecontrol.current_simulation_time+i_fixed_dt, &io_vrt, simVars.benchmark.getExternalForcesUserData);
		simVars.benchmark.getExternalForcesCallback(2, simVars.timecontrol.current_simulation_time+i_fixed_dt, &io_div, simVars.benchmark.getExternalForcesUserData);
	}
}



/*
 * Setup
 */
void SphereAdvection_TS_na_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


std::string SphereAdvection_TS_na_erk::get_help()
{
	std::ostringstream stream;
	stream << " + SphereAdvection_TS_na_erk:" << std::endl;
	stream << "    * 'na_erk'" << std::endl;

	return stream.str();
}


SphereAdvection_TS_na_erk::SphereAdvection_TS_na_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



SphereAdvection_TS_na_erk::~SphereAdvection_TS_na_erk()
{
}

