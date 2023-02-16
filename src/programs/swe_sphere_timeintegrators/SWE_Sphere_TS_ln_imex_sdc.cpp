/*
 * SWE_Sphere_TS_ln_imex_sdc.cpp
 *
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
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
		nNodes(i_simVars.sdc.nNodes),
		nIter(i_simVars.sdc.nIter),
		initialSweepType(i_simVars.sdc.initSweepType),
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
		ts_linear_tendencies_k0(op.sphereDataConfig, i_simVars.sdc.nNodes),
		ts_nonlinear_tendencies_k0(op.sphereDataConfig, i_simVars.sdc.nNodes),
		ts_linear_tendencies_k1(op.sphereDataConfig, i_simVars.sdc.nNodes),
		ts_nonlinear_tendencies_k1(op.sphereDataConfig, i_simVars.sdc.nNodes),

		dt(i_simVars.timecontrol.current_timestep_size)
{
	if (i_simVars.sdc.fileName == "") {
		// No parameter file given, use default SDC settings
		// warning : nNodes default value defined in SimulationVariables.hpp
		//           -> need this to instantiate the tendencies storage
		nIter = 3;
		diagonal = false;
		initialSweepType = "COPY";
		useEndUpdate = false;

		// RADAU-RIGHT nodes, weights quadrature matrix
		tau.setup({nNodes});
		const double tau_default[] = {
			0.15505103, 0.64494897, 1.
		};
		tau = tau_default;

		weights.setup({nNodes});
		const double weights_default[] = {
			0.37640306, 0.51248583, 0.11111111
		};
		weights = weights_default;

		qMat.setup({nNodes, nNodes});
		const double qMat_default[] = {
			0.19681548, -0.06553543, 0.02377097,
			0.39442431,  0.29207341, -0.04154875,
			0.37640306,  0.51248583,  0.11111111
		};
		qMat = qMat_default;

		// BE for implicit sweep
		qMatDeltaI.setup({nNodes, nNodes});
		const double qMatDeltaI_default[] = {
			0.15505103, 0.,         0.,
 			0.15505103, 0.48989795, 0.,
			0.15505103, 0.48989795, 0.        
		};
		qMatDeltaI = qMatDeltaI_default;

		// FE for explicit sweep
		qMatDeltaE.setup({nNodes, nNodes});
		const double qMatDeltaE_default[] = {
			0.,         0.,         0.,
 			0.48989795, 0.,         0.,
			0.48989795, 0.35505103, 0.        
		};
		qMatDeltaE = qMatDeltaE_default;

		// BEpar for initial sweep
		qMatDelta0.setup({nNodes, nNodes});
		const double qMatDelta0_default[] = {
			0.15505103, 0.		  , 0.,
 			0.		  , 0.64494897, 0.,
			0.		  , 0.		  , 1.        
		};
		qMatDelta0 = qMatDelta0_default;
	}
	
	// Initialize irk solver for each nodes
	timestepping_l_irk.resize(nNodes);
	for (int i = 0; i < nNodes; i++)
	{
		timestepping_l_irk[i] = new SWE_Sphere_TS_l_irk(i_simVars, i_op);
		timestepping_l_irk[i]->setup(1, dt*qMatDeltaI(i, i));
	}

	// Print informations ...
	std::cout << "SDC coefficients" << std::endl;
	std::cout << tau << std::endl;
	std::cout << weights << std::endl;
	std::cout << qMat << std::endl;
	std::cout << qMatDeltaE << std::endl;
	std::cout << qMatDeltaI << std::endl;
	std::cout << qMatDelta0 << std::endl;
	std::cout << "nIter :" << nIter << std::endl;
	std::cout << "diagonal : " << diagonal << std::endl;
	std::cout << "useEndUpdate : " << useEndUpdate << std::endl;
	std::cout << "initialSweepType : " << initialSweepType << std::endl;
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
		int iNode
)
{
	timestepping_l_irk[iNode]->solveImplicit(
		rhs.phi, rhs.vrt, rhs.div,
		dt
	);
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
			timestepping_l_irk[i]->update_coefficients(i_fixed_dt*qMatDeltaI(i, i));
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

	// -- compute end-point solution and update step values
	computeEndPoint();

	// -- swap state and time-step
	ts_u0.phi.swapWithConfig(io_phi);
	ts_u0.vrt.swapWithConfig(io_vrt);
	ts_u0.div.swapWithConfig(io_div);
}



void SWE_Sphere_TS_ln_imex_sdc::init_sweep()
{
	// Local convenient references
	const Mat& q = qMat;
	const Mat& qI = qMatDeltaI;
	const Mat& qE = qMatDeltaE;

	if (initialSweepType == "QDELTA")
	{
		// Loop on nodes (can be parallelized if diagonal)

		//SWEET_OMP_PARALLEL_FOR _Pragma("if(diagonal)")
		for (int i = 0; i < nNodes; i++) {

			// Initialize with u0 value
			SWE_VariableVector ts_state(ts_u0);

			if (!diagonal)
			{
				// Add non-linear and linear terms from iteration k (already computed)
				for (int j = 0; j < i; j++) {
					axpy(dt*qE(i, j), ts_nonlinear_tendencies_k0[j], ts_state);
					axpy(dt*qI(i, j), ts_linear_tendencies_k0[j], ts_state);
				}
			}

			// Implicit solve with qI
			solveImplicit(ts_state, dt*qI(i, i), i);

			// Evaluate and store linear term for k
			eval_linear(ts_state, ts_linear_tendencies_k0[i], t0+dt*tau[i]);

			// Evaluate and store non linear term for k
			eval_nonlinear(ts_state, ts_nonlinear_tendencies_k0[i], t0+dt*tau[i]);
		}
	}
	else if (initialSweepType == "COPY")
	{
		// Simply copy to each node as initial guess
		//SWEET_OMP_PARALLEL_FOR
		for (int i = 0; i < nNodes; i++)
		{
			if (i == 0)
			{
				// Evaluate linear and non-linear with initial solution
				eval_linear(ts_u0, ts_linear_tendencies_k0[0], t0);
				eval_nonlinear(ts_u0, ts_nonlinear_tendencies_k0[0], t0);
			}
			else
			{
				ts_linear_tendencies_k0[i] = ts_linear_tendencies_k0[0];
				ts_nonlinear_tendencies_k0[i] = ts_nonlinear_tendencies_k0[0];
			}
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
	//SWEET_OMP_PARALLEL_FOR _Pragma("if(diagonal)")
	for (int i = 0; i < nNodes; i++) {

		// Initialize with u0 value
		SWE_VariableVector ts_state(ts_u0);

		// Add quadrature terms
		for (int j = 0; j < nNodes; j++) {
			axpy(dt*q(i, j), ts_nonlinear_tendencies_k0[j], ts_state);
			axpy(dt*q(i, j), ts_linear_tendencies_k0[j], ts_state);
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

		// Evaluate and store linear term for k+1
		eval_linear(ts_state, ts_linear_tendencies_k1[i], t0+dt*tau[i]);

		// Evaluate and store non linear term for k+1
		eval_nonlinear(ts_state, ts_nonlinear_tendencies_k1[i], t0+dt*tau[i]);
	}

	// Swap k+1 and k values for next iteration (or end-point update)
	ts_nonlinear_tendencies_k0.swap(ts_nonlinear_tendencies_k1);
	ts_linear_tendencies_k0.swap(ts_linear_tendencies_k1);
}

void SWE_Sphere_TS_ln_imex_sdc::computeEndPoint()
{
	if (useEndUpdate)
	{
		/*
		 * Use quadrature
		 */
		const Vec& w = weights;
		
		// Compute collocation update
		SWE_VariableVector ts_state(ts_u0);

		for (int j = 0; j < nNodes; j++) {
			axpy(dt*w(j), ts_nonlinear_tendencies_k0[j], ts_state);
			axpy(dt*w(j), ts_linear_tendencies_k0[j], ts_state);
		}

		// Time-step update using last state value
		ts_u0.phi = ts_state.phi;
		ts_u0.vrt = ts_state.vrt;
		ts_u0.div = ts_state.div;
		return;
	}

	SWEETError("TODO");
#if 0
	/*
	 * Use final point in time
	 */
	ts_u0 = ts_nonlinear_tendencies_k0[nNodes-1];
#endif
}

/*
 * Setup
 */
void SWE_Sphere_TS_ln_imex_sdc::setup()
{
	for (int i = 0; i < nNodes; i++)
	{
		timestepping_l_irk[i]->setup(
			1,
			dt*qMatDeltaI(i, i)
		);
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
	}
	
}

