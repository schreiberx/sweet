/*
 * Burgers_Plane.hpp
 *
 *  Created on: 15 May 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */
#ifndef SRC_PROGRAMS_BURGERSPLANE_HPP_
#define SRC_PROGRAMS_BURGERSPLANE_HPP_


#include <complex>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>


#if SWEET_MPI
#	include <mpi.h>
#endif

/**
 * This class implements the REXI (rational approximation of exponential integrator) solver for the SWE,
 * see High-order time-parallel approximation of evolution operators, T. Haut et al.
 *
 * We split this file into the header and cpp file.
 *
 * This allows using a OpenMP parallelization only in the RexiSWE class to check this degree of
 * parallelism.
 */
class Burgers_Plane
{

	/**
	 * IMEX time stepping for the coarse timestepping method
	 *
	 * The IMEX RK schemes combine an implicit and explicit time discretization.
	 * The diffusive, stiff term gets treated implicitly and the convective, non-
	 * stiff term gets treated explicitly. This results in the following system,
	 * which is solved by this routine:
	 * (I-\nu\Delta) u(t+\tau) = u(t) - \tau (u(t)*\nabla(u(t))
	 * for u(t+\tau).
	 */
public:
	bool run_timestep_imex(
		PlaneData &io_u,
		PlaneData &io_v,

		double& o_dt,			///< return time step size for the computed time step
		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_parameters,
		PlaneDataConfig* planeDataConfig,

		double i_max_simulation_time = std::numeric_limits<double>::infinity()	///< limit the maximum simulation time
	);


	/**
	 * Solve U_t = L U via implicit solver:
	 *
	 */
public:
	bool run_timestep_implicit_ts(
		PlaneData &io_u,
		PlaneData &io_v,

		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_parameters
	);


	/**
	 * Solve U_t = L U via explicit solver:
	 *
	 */
public:
	bool run_timestep_explicit_ts(
			PlaneData &io_u,
			PlaneData &io_v,

			double& o_dt,			///< return time step size for the computed time step
			double i_timestep_size,	///< timestep size

			PlaneOperators &op,
			const SimulationVariables &i_parameters,
			PlaneDataConfig* planeDataConfig
	);


	/**
	 * Solve  viscous Burgers' equation with semi-Lagrangian
	 *
	 */
public:
	bool run_timestep_sl(
		PlaneData &io_u,  ///< Current and past fields
		PlaneData &io_v,
		PlaneData &io_u_prev,
		PlaneData &io_v_prev,

		ScalarDataArray &i_posx_a, //Arrival point positions in x and y (this is basically the grid)
		ScalarDataArray &i_posy_a,

		double& o_dt,			///< return time step size for the computed time step
		double i_timestep_size,	///< timestep size

		const SimulationVariables &i_simVars, ///< Parameters for simulation
		PlaneDataConfig* planeDataConfig,

		PlaneOperators &op,     ///< Operator class
		PlaneDataSampler &sampler2D, ///< Interpolation class
		SemiLagrangian &semiLagrangian,  ///< Semi-Lag class
		double* i_stag_displacement,
		double* i_stag_u,
		double* i_stag_v
	);

	/**
	 * This method computes the analytical solution based on the given initial values.
	 */
public:
	void run_timestep_direct_solution(
			PlaneData &io_h,
			PlaneData &io_u,
			PlaneData &io_v,

			double i_timestep_size,	///< timestep size

			PlaneOperators &op,
			const SimulationVariables &i_parameters
	);

};

#endif /* SRC_PROGRAMS_BURGERSPLANE_HPP_ */
