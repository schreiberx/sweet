/*
 * SWE_Sphere_TS_ln_imex_sdc.cpp
 *
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_ln_imex_sdc.hpp"


bool SWE_Sphere_TS_ln_imex_sdc::implements_timestepping_method(const std::string &i_timestepping_method
									)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;

	if (i_timestepping_method == "ln_imex_sdc")
		return true;

	return false;
}



void SWE_Sphere_TS_ln_imex_sdc::setup_auto()
{
	setup(timestepping_order, timestepping_order2, 0);
}



std::string SWE_Sphere_TS_ln_imex_sdc::string_id()
{
	return "ln_imex_sdc";
}


void SWE_Sphere_TS_ln_imex_sdc::run_timestep(
		SphereData_Spectral &io_phi,
		SphereData_Spectral &io_vrt,
		SphereData_Spectral &io_div,
		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	/*
	Solves the following problem
	
	du/dt = L(u) + NL(u)

	with u = [io_phi, io_vrt, io_div], L the linear term
	and NL the non linear term.
	u_in is the initial value before any method call,
	and u_out is the final value after.
	*/

	// first order explicit for non-linear
	// -- u_out = u_in + dt*NL(u_in,t)
	timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
			io_phi, io_vrt, io_div,
			i_fixed_dt,
			i_simulation_timestamp
		);

	// first order IRK for linear
	// -- u_out - dt*L(u_out, t) = u_in ... actually, not really ...
	timestepping_l_irk.run_timestep(
			io_phi, io_vrt, io_div,
			i_fixed_dt,
			i_simulation_timestamp
		);

	// TODO : Implement IMEX SDC instead !
}

void SWE_Sphere_TS_ln_imex_sdc::evalLinearTerms(
		SphereData_Spectral &phi_pert,	///< prognostic variables
		SphereData_Spectral &vort,	    ///< prognostic variables
		SphereData_Spectral &div,	    ///< prognostic variables
		SphereData_Spectral &phi_pert_L,	///< evaluation
		SphereData_Spectral &vort_L,	    ///< evaluation
		SphereData_Spectral &div_L,	        ///< evaluation
		double simulation_timestamp
) {
	timestepping_l_erk_n_erk.euler_timestep_update_linear(
		phi_pert,
		vort,
		div,
		// Variables where to store the linear terms evaluation
		phi_pert_L,
		vort_L,
		div_L,
		// Timestamp of the evaluation
		simulation_timestamp
	);
}

void SWE_Sphere_TS_ln_imex_sdc::evalNonLinearTerms(
		SphereData_Spectral &phi_pert,	///< prognostic variables
		SphereData_Spectral &vort,	    ///< prognostic variables
		SphereData_Spectral &div,	    ///< prognostic variables
		SphereData_Spectral &phi_pert_N,	///< evaluation
		SphereData_Spectral &vort_N,	    ///< evaluation
		SphereData_Spectral &div_N,	    ///< evaluation
		double simulation_timestamp
) {
	timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
		phi_pert,
		vort,
		div,
		// Variables where to store the non linear terms evaluation
		phi_pert_N,
		vort_N,
		div_N,
		// Timestamp of the evaluation
		simulation_timestamp
	);
}

void SWE_Sphere_TS_ln_imex_sdc::solveImplicit(
		SphereData_Spectral &rhs_phi,	///< rhs variables
		SphereData_Spectral &rhs_vrt,	///< rhs variables
		SphereData_Spectral &rhs_div,	///< rhs variables

		double dt
) {
	timestepping_l_irk.solveImplicit(
		rhs_phi,
		rhs_vrt,
		rhs_div,
		dt
	);
}

/*
 * Setup
 */
void SWE_Sphere_TS_ln_imex_sdc::setup(
		int i_order,	///< order of RK time stepping method
		int i_order2,	///< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	if (i_order2 < 0)
		i_order2 = i_order;

	if (i_order != i_order2)
		SWEETError("Orders of 1st and 2nd one must match");

	version_id = i_version_id;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;
	timestep_size = simVars.timecontrol.current_timestep_size;

	if (timestepping_order == 1)
	{
		timestepping_l_irk.setup(
			1,
			timestep_size
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_l_irk.setup(
					2,
					timestep_size*0.5,
					simVars.disc.timestepping_crank_nicolson_filter,
					false
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_irk.setup(
					2,
					timestep_size,
					simVars.disc.timestepping_crank_nicolson_filter,
					false
			);
		}
		else
		{
			SWEETError("Invalid version");
		}
	}
	else
	{
		SWEETError("Invalid timestepping order");
	}


	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_l_erk_n_erk.setup(1, 1);
}



SWE_Sphere_TS_ln_imex_sdc::SWE_Sphere_TS_ln_imex_sdc(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_l_irk(simVars, op),
		timestepping_l_erk_n_erk(simVars, op),
		version_id(0),
		timestepping_order(-1)
{

}



SWE_Sphere_TS_ln_imex_sdc::~SWE_Sphere_TS_ln_imex_sdc()
{
}

