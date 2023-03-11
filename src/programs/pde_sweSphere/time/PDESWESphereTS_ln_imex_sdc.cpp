/*
 * PDESWESphereTS_ln_imex_sdc.cpp
 *
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_ln_imex_sdc.hpp"




bool PDESWESphereTS_ln_imex_sdc::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup(io_ops, shackPDESWETimeDisc->timestepping_order, shackPDESWETimeDisc->timestepping_order2, 0);
}


/*
 * Setup
 */
bool PDESWESphereTS_ln_imex_sdc::setup(
		sweet::SphereOperators *io_ops,
		int i_order,	///< order of RK time stepping method
		int i_order2,	///< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	ops = io_ops;
	timestepping_order = i_order;
	timestepping_order2 = i_order2;
	timestep_size = shackTimestepControl->current_timestep_size;


	if (timestepping_order2 < 0)
		timestepping_order2 = timestepping_order;

	if (timestepping_order != timestepping_order2)
		SWEETError("Orders of 1st and 2nd one must match");

#if 0
	/*
	 * TODO
	 */

	// SDC main parameters
	nNodes = i_shackDict.sdc.nNodes;
	nIter(i_shackDict.sdc.nIter),
	initialSweepType(i_shackDict.sdc.initSweepType),
	diagonal(i_shackDict.sdc.diagonal),
	useEndUpdate(i_shackDict.sdc.useEndUpdate),

	// Nodes and weights
	tau(i_shackDict.sdc.nodes),
	weights(i_shackDict.sdc.weights),

	// Quadrature matrices
	qMat(i_shackDict.sdc.qMatrix),
	qMatDeltaI(i_shackDict.sdc.qDeltaI),
	qDeltaE(i_shackDict.sdc.qDeltaE),
	qMatDelta0(i_shackDict.sdc.qDelta0),

	// Storage and reference container (for u0)
	ts_linear_tendencies_k0(ops->sphereDataConfig, i_shackDict.sdc.nNodes),
	ts_nonlinear_tendencies_k0(ops->sphereDataConfig, i_shackDict.sdc.nNodes),
	ts_linear_tendencies_k1(ops->sphereDataConfig, i_shackDict.sdc.nNodes),
	ts_nonlinear_tendencies_k1(ops->sphereDataConfig, i_shackDict.sdc.nNodes),
	ts_tmp_state(ops->sphereDataConfig)
#endif


	timestepping_l_irk.setup(
		ops,
		1,
		timestep_size
	);


	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_l_erk_n_erk.setup_main(ops, 1, 1);

	return true;
}



bool PDESWESphereTS_ln_imex_sdc::implementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;
	if (i_timestepping_method == "ln_imex_sdc")
		return true;

	return false;
}




std::string PDESWESphereTS_ln_imex_sdc::getIDString()
{
	return "ln_imex_sdc";
}


void PDESWESphereTS_ln_imex_sdc::runTimestep(
		sweet::SphereData_Spectral &io_phi,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,
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

	// Default behavior (no parameter set ...)
	if (nNodes == 0) {
		return;
	}

	// -- set-up step values container
	ts_u0.phi.swapWithConfig(io_phi);
	ts_u0.vrt.swapWithConfig(io_vrt);
	ts_u0.div.swapWithConfig(io_div);

	//u0.setPointer(io_phi, io_vrt, io_div);
	t0 = i_simulation_timestamp;
	dt = i_fixed_dt;

	// -- initialize nodes values and state
	init_sweep();

	// -- perform sweeps
	for (int k = 0; k < nIter; k++){
		sweep(k);
	}

	// -- compute end-point solution and update step values
	computeEndPoint();

	ts_u0.phi.swapWithConfig(io_phi);
	ts_u0.vrt.swapWithConfig(io_vrt);
	ts_u0.div.swapWithConfig(io_div);
}

void PDESWESphereTS_ln_imex_sdc::eval_linear(const SWE_VariableVector& u, SWE_VariableVector& eval, double t) {
	timestepping_l_erk_n_erk.euler_timestep_update_linear(
		u.phi, u.vrt, u.div,
		eval.phi, eval.vrt, eval.div,
		t  // time-stamp
	);
}
void PDESWESphereTS_ln_imex_sdc::eval_nonlinear(const SWE_VariableVector& u, SWE_VariableVector& eval, double t) {
	timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
		u.phi, u.vrt, u.div,
		eval.phi, eval.vrt, eval.div,
		t  // time-stamp
	);
}

void PDESWESphereTS_ln_imex_sdc::solveImplicit(
		SWE_VariableVector& rhs,
		double dt
)
{
	timestepping_l_irk.solveImplicit(
		rhs.phi, rhs.vrt, rhs.div,
		dt
	);
}

void PDESWESphereTS_ln_imex_sdc::axpy(
		double a,
		const SWE_VariableVector& i_x,
		SWE_VariableVector& io_y
)
{
	io_y.phi += a*i_x.phi;
	io_y.vrt += a*i_x.vrt;
	io_y.div += a*i_x.div;
}

void PDESWESphereTS_ln_imex_sdc::init_sweep() {
	
	if (initialSweepType == "QDELTA")
	{
		// Uses QDelta matrix to initialize node values

		// Local convenient references
		//Mat& q = qMat;
		Mat& qI = qMatDeltaI;
		Mat& qE = qMatDeltaI;
		Mat& q0 = qMatDelta0;

		// Loop on nodes (can be parallelized if diagonal)
		for (int i = 0; i < nNodes; i++) {

			// Initialize with u0 value
			ts_tmp_state = ts_u0;

			if (diagonal) {
				// qMatDelta0 is diagonal-only
				// Hence, Here, we only iterate over the diagonal
				// Implicit solve with q0
				solveImplicit(ts_tmp_state, dt*q0(i, i));
			} else {
				// Add non-linear and linear terms from iteration k (already computed)
				for (int j = 0; j < i; j++) {
					axpy(dt*qE(i, j), ts_nonlinear_tendencies_k0[j], ts_tmp_state);
					axpy(dt*qI(i, j), ts_linear_tendencies_k0[j], ts_tmp_state);
				}
				// Implicit solve with qI
				solveImplicit(ts_tmp_state, dt*qI(i, i));
			}

			// Evaluate and store linear term for k
			eval_linear(ts_tmp_state, ts_linear_tendencies_k0[i], t0+dt*tau(i));

			// Evaluate and store non linear term for k
			eval_nonlinear(ts_tmp_state, ts_nonlinear_tendencies_k0[i], t0+dt*tau(i));
		}

	} else if (initialSweepType == "COPY") {
		// Evaluate linear and non-linear with initial solution
		eval_linear(ts_u0, ts_linear_tendencies_k0[0], t0);
		eval_nonlinear(ts_u0, ts_nonlinear_tendencies_k0[0], t0);

		// Simply copy to each node as initial guess
		for (int i = 1; i < nNodes; i++) {
			ts_linear_tendencies_k0[i] = ts_linear_tendencies_k0[0];
			ts_nonlinear_tendencies_k0[i] = ts_nonlinear_tendencies_k0[0];
		}
	} else {
		SWEETError(initialSweepType + " initial sweep type not implemented");
	}

}

void PDESWESphereTS_ln_imex_sdc::sweep(size_t k) {

	// Local convenient references
	const Mat& q = qMat;
	const Mat& qDeltaImplicit = qMatDeltaI;
	const Mat& qDeltaExplicit = qMatDeltaI;

	// Loop on nodes (can be parallelized if diagonal)
	for (int i = 0; i < nNodes; i++) {
		
		// Initialize with u0 value
		ts_tmp_state = ts_u0;
		
		// Add quadrature terms
		for (int j = 0; j < nNodes; j++) {
			//double a = q(i, j);
			axpy(dt*q(i, j), ts_nonlinear_tendencies_k0[j], ts_tmp_state);
			axpy(dt*q(i, j), ts_linear_tendencies_k0[j], ts_tmp_state);
		}

		if (!diagonal) {
			// Add non-linear and linear terms from iteration k+1
			for (int j = 0; j < i; j++) {
				axpy(dt*qDeltaExplicit(i, j), ts_nonlinear_tendencies_k1[j], ts_tmp_state);
				axpy(dt*qDeltaImplicit(i, j), ts_linear_tendencies_k1[j], ts_tmp_state);
			}

			// Substract non-linear and linear terms from iteration k
			for (int j = 0; j < i; j++) {
				axpy(-dt*qDeltaExplicit(i, j), ts_nonlinear_tendencies_k0[j], ts_tmp_state);
				axpy(-dt*qDeltaImplicit(i, j), ts_linear_tendencies_k0[j], ts_tmp_state);
			}
		}

		// Last linear term from iteration k
		axpy(-dt*qDeltaImplicit(i, i), ts_linear_tendencies_k0[i], ts_tmp_state);
		
		// Implicit solve
		solveImplicit(ts_tmp_state, dt*qDeltaImplicit(i,i));

		// Evaluate and store linear term for k+1
		eval_linear(ts_tmp_state, ts_linear_tendencies_k1[i], t0+dt*tau(i));

		// Evaluate and store non linear term for k+1
		eval_nonlinear(ts_tmp_state, ts_nonlinear_tendencies_k1[i], t0+dt*tau(i));
	}

	// Swap k+1 and k values for next iteration (or end-point update)

	ts_nonlinear_tendencies_k0.swap(ts_nonlinear_tendencies_k1);
	ts_linear_tendencies_k0.swap(ts_linear_tendencies_k1);
}

void PDESWESphereTS_ln_imex_sdc::computeEndPoint()
{
	if (useEndUpdate) {
		// Compute collocation update
		const Vec& w = weights;
		ts_tmp_state = ts_u0;
		for (int j = 0; j < nNodes; j++) {
			axpy(dt*w(j), ts_nonlinear_tendencies_k0[j], ts_tmp_state);
			axpy(dt*w(j), ts_linear_tendencies_k0[j], ts_tmp_state);
		}
	}

	// Time-step update using last state value
	ts_u0.phi = ts_tmp_state.phi;
	ts_u0.vrt = ts_tmp_state.vrt;
	ts_u0.div = ts_tmp_state.div;
}



PDESWESphereTS_ln_imex_sdc::PDESWESphereTS_ln_imex_sdc()
{
}



PDESWESphereTS_ln_imex_sdc::~PDESWESphereTS_ln_imex_sdc()
{
}

