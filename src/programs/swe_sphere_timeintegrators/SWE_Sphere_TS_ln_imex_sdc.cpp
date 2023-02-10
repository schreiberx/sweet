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
	Performe one step of IMEX SDC for
	
	du/dt = L(u) + NL(u)

	with u = [io_phi, io_vrt, io_div], L the linear term
	and NL the non linear term.

	It uses an implicit sweep for the L term,
	and an explicit sweep for the NL term.
	*/

	// -- set-up step values container
	u0.setRef(io_phi, io_vrt, io_div);
	t0 = i_simulation_timestamp;
	dt = i_fixed_dt;

	// -- initialize nodes values and state
	initSweep();

	// -- perform sweeps
	for (size_t k = 0; k < nIter; k++){
		sweep(k);
	}

	// -- compute end-point solution and update step values
	prolongate();
}

void SWE_Sphere_TS_ln_imex_sdc::evalLinearTerms(const SWE_Variables& u, SWE_Variables& eval, double t) {
	timestepping_l_erk_n_erk.euler_timestep_update_linear(
		u.phi, u.vort, u.div,
		eval.phi, eval.vort, eval.div,
		t  // time-stamp
	);
}

void SWE_Sphere_TS_ln_imex_sdc::evalNonLinearTerms(const SWE_Variables& u, SWE_Variables& eval, double t) {
	timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
		u.phi, u.vort, u.div,
		eval.phi, eval.vort, eval.div,
		t  // time-stamp
	);
}

void SWE_Sphere_TS_ln_imex_sdc::solveImplicit(SWE_Variables& rhs, double dt) {
	timestepping_l_irk.solveImplicit(
		rhs.phi, rhs.vort, rhs.div,
		dt
	);
}

void SWE_Sphere_TS_ln_imex_sdc::axpy(double a, const SWE_Variables& x, SWE_Variables& y) {
	// TODO : optimize this if possible (check with Martin)
	tmp = x.phi;
	tmp *= a;
	y.phi += tmp;

	tmp = x.vort;
	tmp *= a;
	y.vort += tmp;

	tmp = x.div;
	tmp *= a;
	y.div += tmp;
}

void SWE_Sphere_TS_ln_imex_sdc::initSweep() {
	// Initialize state with step values
	state.fillWith(u0);

	// Evaluate linear and non-linear with initial solution, and copy to each node
	evalLinearTerms(state, lTerms.getK(0), t0);
	evalNonLinearTerms(state, nTerms.getK(0), t0);
	for (size_t i = 1; i < nNodes; i++) {
		lTerms.getK(i).fillWith(lTerms.getK(0));
		nTerms.getK(i).fillWith(nTerms.getK(0));
	}
}

void SWE_Sphere_TS_ln_imex_sdc::sweep(size_t k) {

	// Local convenient references
	const Mat& q = qMat;
	const Mat& qI = qDeltaI;
	const Mat& qE = qDeltaI;

	// Loop on nodes
	for (size_t i = 0; i < nNodes; i++) {
		
		// Initialize with u0 value
		state.fillWith(u0);
		
		// Add quadrature terms
		for (size_t j = 0; j < nNodes; j++) {
			axpy(dt*q[i][j], nTerms.getK(j), state);
			axpy(dt*q[i][j], lTerms.getK(j), state);
		}

		// // Add non-linear and linear terms from iteration k+1
		for (size_t j = 0; j < i; j++) {
			axpy(dt*qE[i][j], nTerms.getK1(j), state);
			axpy(dt*qI[i][j], lTerms.getK1(j), state);
		}

		// Substract non-linear and linear terms from iteration k
		for (size_t j = 0; j < i; j++) {
			axpy(-dt*qE[i][j], nTerms.getK(j), state);
			axpy(-dt*qI[i][j], lTerms.getK(j), state);
		}
		axpy(-dt*qI[i][i], lTerms.getK(i), state);
		
		// Implicit solve
		// solveImplicit(state, dt*qI[i][i]);

		// Evaluate and store linear term for k+1
		evalLinearTerms(state, lTerms.getK1(i), t0+dt*tau[i]);

		// Evaluate and store non linear term for k+1
		evalNonLinearTerms(state, nTerms.getK1(i), t0+dt*tau[i]);
	}

	// Swap k+1 and k values for next iteration
	nTerms.swapIterate();
	lTerms.swapIterate();
}

void SWE_Sphere_TS_ln_imex_sdc::prolongate() {
	*(u0.phi) = state.phi;
	*(u0.vort) = state.vort;
	*(u0.div) = state.div;
	// TODO : implement quadrature prolongation
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
		timestepping_order(-1),
		lTerms(op.sphereDataConfig),
		nTerms(op.sphereDataConfig),
		state(op.sphereDataConfig),
		tmp(op.sphereDataConfig),
		u0()
{

}



SWE_Sphere_TS_ln_imex_sdc::~SWE_Sphere_TS_ln_imex_sdc()
{
}

