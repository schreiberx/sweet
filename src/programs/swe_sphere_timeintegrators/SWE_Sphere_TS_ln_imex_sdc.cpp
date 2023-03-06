/*
 * SWE_Sphere_TS_ln_imex_sdc.cpp
 *
 *      Author: Thibaut LUNET <thibaut.lunet@tuhh.de>
 */

#include "SWE_Sphere_TS_ln_imex_sdc.hpp"

#include <sweet/TimeStepSizeChanged.hpp>


SWE_Sphere_TS_ln_imex_sdc::SWE_Sphere_TS_ln_imex_sdc(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_l_erk_n_erk(simVars, op),
		timestepping_order(-1),

		// SDC main parameters
		nNodes(i_simVars.sdc.nodes.size()),
		nIter(i_simVars.sdc.nIter),
		initialSweepType(i_simVars.sdc.initialSweepType),
		diagonal(i_simVars.sdc.diagonal),
		useEndUpdate(i_simVars.sdc.useEndUpdate),

		// Nodes and weights
		tau(i_simVars.sdc.nodes),
		weights(i_simVars.sdc.weights),

		// Quadrature matrices
		qMat(i_simVars.sdc.qMatrix),
		qMatDeltaI(i_simVars.sdc.qDeltaI),
		qMatDeltaE(i_simVars.sdc.qDeltaE),
		qMatDelta0(i_simVars.sdc.qDelta0),

		// Storage and reference container (for u0)
		ts_linear_tendencies_k0(op.sphereDataConfig, i_simVars.sdc.nodes.size()),
		ts_nonlinear_tendencies_k0(op.sphereDataConfig, i_simVars.sdc.nodes.size()),
		ts_linear_tendencies_k1(op.sphereDataConfig, i_simVars.sdc.nodes.size()),
		ts_nonlinear_tendencies_k1(op.sphereDataConfig, i_simVars.sdc.nodes.size()),

		dt(i_simVars.timecontrol.current_timestep_size),
		idString(i_simVars.sdc.idString)
{
	
	// Initialize irk solver for each nodes
	timestepping_l_irk.resize(nNodes);
	if (initialSweepType == "QDELTA") {
		timestepping_l_irk_init.resize(nNodes);
	}

#if SWEET_PARALLEL_SDC_OMP_MODEL
	SWEET_OMP_PARALLEL_FOR _Pragma("num_threads(nNodes)")
#endif
	for (int i = 0; i < nNodes; i++)
	{
		// Initialize LHS coefficients for sweeps
		timestepping_l_irk[i] = new SWE_Sphere_TS_l_irk(i_simVars, i_op);
		timestepping_l_irk[i]->setup(1, dt*qMatDeltaI(i, i));
		
		if (initialSweepType == "QDELTA") {
			// Initialize LHS coefficients for initial sweep with QDELTA
			timestepping_l_irk_init[i] = new SWE_Sphere_TS_l_irk(i_simVars, i_op);
			timestepping_l_irk_init[i]->setup(1, dt*qMatDelta0(i, i));
		}
	}

}


bool SWE_Sphere_TS_ln_imex_sdc::implements_timestepping_method(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;
	if (i_timestepping_method == "ln_imex_sdc")
		return true;
	return false;
}



void SWE_Sphere_TS_ln_imex_sdc::setup_auto()
{
	setup();
}



std::string SWE_Sphere_TS_ln_imex_sdc::string_id()
{
	return "ln_imex_sdc";
}

void SWE_Sphere_TS_ln_imex_sdc::eval_linear(const SWE_VariableVector& i_u, SWE_VariableVector& o_eval, double t) {
	timestepping_l_erk_n_erk.euler_timestep_update_linear(
		i_u.phi, i_u.vrt, i_u.div,
		o_eval.phi, o_eval.vrt, o_eval.div,
		t  // time-stamp
	);
}
void SWE_Sphere_TS_ln_imex_sdc::eval_nonlinear(const SWE_VariableVector& i_u, SWE_VariableVector& o_eval, double t) {
	timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
		i_u.phi, i_u.vrt, i_u.div,
		o_eval.phi, o_eval.vrt, o_eval.div,
		t  // time-stamp
	);
}

void SWE_Sphere_TS_ln_imex_sdc::solveImplicit(
		SWE_VariableVector& rhs,
		double dt,
		int iNode,
		bool initSweep
)
{
	if (initSweep) 
	{
		timestepping_l_irk_init[iNode]->solveImplicit(
			rhs.phi, rhs.vrt, rhs.div,
			dt
		);
	} else {
		timestepping_l_irk[iNode]->solveImplicit(
			rhs.phi, rhs.vrt, rhs.div,
			dt
		);
	}
}

void SWE_Sphere_TS_ln_imex_sdc::axpy(
		double a,
		const SWE_VariableVector& i_x,
		SWE_VariableVector& io_y
)
{
	io_y.phi += a*i_x.phi;
	io_y.vrt += a*i_x.vrt;
	io_y.div += a*i_x.div;
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
	Perform one step of IMEX SDC for
	
	du/dt = L(u) + NL(u)

	with u = [io_phi, io_vrt, io_div], L the linear term
	and NL the non linear term.

	It uses an implicit sweep for the L term,
	and an explicit sweep for the NL term.
	*/

	// -- set-up step values container
	ts_u0.phi.swapWithConfig(io_phi);
	ts_u0.vrt.swapWithConfig(io_vrt);
	ts_u0.div.swapWithConfig(io_div);

	if (TimeStepSizeChanged::is_changed(dt, i_fixed_dt, true)){
		std::cout << "SDC: UPDATING LHS COEFFICIENTS" << std::endl;
		for (int i = 0; i < nNodes; i++)
		{
			// LHS coefficients for sweeps
			timestepping_l_irk[i]->update_coefficients(i_fixed_dt*qMatDeltaI(i, i));

			if (initialSweepType == "QDELTA") {
				// LHS coefficients for initial sweep
				timestepping_l_irk_init[i]->update_coefficients(i_fixed_dt*qMatDelta0(i, i));
			}
		}
	}
	t0 = i_simulation_timestamp;
	dt = i_fixed_dt;
	

	// -- initialize nodes values and state
	init_sweep();

	// -- perform sweeps
	for (int k = 0; k < nIter; k++){
		sweep(k);
	}

	// -- swap state and time-step
	ts_u0.phi.swapWithConfig(io_phi);
	ts_u0.vrt.swapWithConfig(io_vrt);
	ts_u0.div.swapWithConfig(io_div);
}



void SWE_Sphere_TS_ln_imex_sdc::init_sweep()
{
	// Local convenient references
	const Mat& q = qMat;
	const Mat& q0 = qMatDelta0;
	const Mat& qE = qMatDeltaE;

	if (initialSweepType == "QDELTA")
	{

		// Loop on nodes (can be parallelized if diagonal)

#if SWEET_PARALLEL_SDC_OMP_MODEL
		#pragma omp parallel num_threads(nNodes) if(diagonal)
		#pragma omp for schedule(static,1)
#endif
		for (int i = 0; i < nNodes; i++) {

			// Initialize with u0 value
			SWE_VariableVector ts_state(ts_u0);

			if (!diagonal)
			{
				// Add non-linear and linear terms from iteration k (already computed)
				for (int j = 0; j < i; j++) {
					axpy(dt*qE(i, j), ts_nonlinear_tendencies_k0[j], ts_state);
					axpy(dt*q0(i, j), ts_linear_tendencies_k0[j], ts_state);
				}
			}

			// Implicit solve with qI
			solveImplicit(ts_state, dt*q0(i, i), i, true);

			double tNode = t0+dt*tau(i);

			// Evaluate and store linear term for k
			eval_linear(ts_state, ts_linear_tendencies_k0[i], tNode);

			// Evaluate and store non linear term for k
			eval_nonlinear(ts_state, ts_nonlinear_tendencies_k0[i], tNode);
		}
	}
	else if (initialSweepType == "COPY")
	{
		// Evaluate linear and non-linear with initial solution
		eval_linear(ts_u0, ts_linear_tendencies_k0[0], t0);
		eval_nonlinear(ts_u0, ts_nonlinear_tendencies_k0[0], t0);

		// Simply copy to each node as initial guess
#if SWEET_PARALLEL_SDC_OMP_MODEL
		#pragma omp parallel num_threads(nNodes)
		#pragma omp for schedule(static,1)
#endif
		for (int i = 0; i < nNodes; i++)
		{
			// Include first node just to simulate same parallelization as for parallel updates
			if (i == 0)
				continue;

			ts_linear_tendencies_k0[i] = ts_linear_tendencies_k0[0];
			ts_nonlinear_tendencies_k0[i] = ts_nonlinear_tendencies_k0[0];
		}
	}
	else
	{
		SWEETError(initialSweepType + " initial sweep type not implemented");
	}

}


void SWE_Sphere_TS_ln_imex_sdc::sweep(
		size_t k	/// iteration number
) {

	// Local convenient references
	const Mat& q = qMat;
	const Mat& qI = qMatDeltaI;
	const Mat& qE = qMatDeltaE;

	// Loop on nodes (can be parallelized if diagonal)
#if SWEET_PARALLEL_SDC_OMP_MODEL
	#pragma omp parallel num_threads(nNodes) if(diagonal)
	#pragma omp for schedule(static,1)
#endif
	for (int i = 0; i < nNodes; i++) {

		// Initialize with u0 value
		SWE_VariableVector ts_state(ts_u0);

		// Add quadrature terms
		for (size_t j = 0; j < nNodes; j++) {
			double a = dt*q(i, j);
			axpy(a, ts_nonlinear_tendencies_k0[j], ts_state);
			axpy(a, ts_linear_tendencies_k0[j], ts_state);
		}

		if (!diagonal) {

			// Add non-linear and linear terms from iteration k+1
			for (int j = 0; j < i; j++) {
				axpy(dt*qE(i, j), ts_nonlinear_tendencies_k1[j], ts_state);
				axpy(dt*qI(i, j), ts_linear_tendencies_k1[j], ts_state);
			}

			// Substract non-linear and linear terms from iteration k
			for (int j = 0; j < i; j++) {
				axpy(-dt*qE(i, j), ts_nonlinear_tendencies_k0[j], ts_state);
				axpy(-dt*qI(i, j), ts_linear_tendencies_k0[j], ts_state);
			}
		}

		// Last linear term from iteration k
		axpy(-dt*qI(i, i), ts_linear_tendencies_k0[i], ts_state);

		// Implicit solve
		solveImplicit(ts_state, dt*qI(i,i), i);

		if (!useEndUpdate && (i == nNodes-1) && (k == nIter-1)) {
			// Time-step update using last state value
			ts_u0.phi = ts_state.phi;
			ts_u0.vrt = ts_state.vrt;
			ts_u0.div = ts_state.div;
		} 
		else
		{

			double tNode = t0+dt*tau(i);

			// Evaluate and store linear term for k+1
			eval_linear(ts_state, ts_linear_tendencies_k1[i], tNode);

			// Evaluate and store non linear term for k+1
			eval_nonlinear(ts_state, ts_nonlinear_tendencies_k1[i], tNode);
		}
	}

	if (!useEndUpdate && (k == nIter-1))
		return;

	// Swap k+1 and k values for next iteration (or end-point update)
	ts_nonlinear_tendencies_k0.swap(ts_nonlinear_tendencies_k1);
	ts_linear_tendencies_k0.swap(ts_linear_tendencies_k1);

	// Compute end point for last sweep
	if (k+1 == nIter) 
	{
		if (useEndUpdate) {

			/*
			* Use quadrature
			*/
			const Vec& w = weights;
			
			// Compute collocation update
			SWE_VariableVector ts_state(ts_u0);

			for (int j = 0; j < nNodes; j++) {
				double a = dt*w(j);
				axpy(a, ts_nonlinear_tendencies_k0[j], ts_state);
				axpy(a, ts_linear_tendencies_k0[j], ts_state);
			}
			// Time-step update using last state value
			ts_u0.phi = ts_state.phi;
			ts_u0.vrt = ts_state.vrt;
			ts_u0.div = ts_state.div;
		}
	}
}

/*
 * Setup
 */
void SWE_Sphere_TS_ln_imex_sdc::setup()
{
	for (int i = 0; i < nNodes; i++)
	{
		timestepping_l_irk[i]->setup(1, dt*qMatDeltaI(i, i));
		if (initialSweepType == "QDELTA") {
			timestepping_l_irk_init[i]->setup(1, dt*qMatDelta0(i, i));
		}
		
	}

	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_l_erk_n_erk.setup(1, 1);
}


SWE_Sphere_TS_ln_imex_sdc::~SWE_Sphere_TS_ln_imex_sdc()
{
	for (int i = 0; i < nNodes; i++)
	{
		delete timestepping_l_irk[i];
		if (initialSweepType == "QDELTA") {
			delete timestepping_l_irk_init[i];
		}
	}
	
}

