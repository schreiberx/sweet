/*
 * SWE_Plane_TS_l_irk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_IRK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_IRK_HPP_

#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include "SWE_Plane_TS_interface.hpp"



class SWE_Plane_TS_l_irk	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;
	PlaneOperatorsComplex opComplex;

	int timestepping_order;

public:
	SWE_Plane_TS_l_irk(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			PlaneData &io_h,	///< prognostic variables
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	);


	/**
	 * Solve complex-valued Helmholtz problem with a spectral solver,
	 * values are given in Spectral space
	 *
	 * (kappa - gh*D2) X = B
	 *
	 * This is here to compare speedups and such a solver cannot be applied in general,
	 * e.g. with special general boundary values
	 */
public:
	void helmholtz_spectral_solver_spec(
			std::complex<double> i_kappa,
			double i_gh0,
			const PlaneDataComplex &i_rhs,
			PlaneDataComplex &io_x,
			int i_thread_id = 0
	)
	{
		PlaneDataComplex lhs = (-i_gh0*(opComplex.diff2_c_x + opComplex.diff2_c_y)).spectral_addScalarAll(i_kappa);

		io_x = i_rhs.spectral_div_element_wise(lhs);
	}

        PlaneOperatorsComplex& get_plane_operators_complex() {return opComplex;}

	virtual ~SWE_Plane_TS_l_irk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
