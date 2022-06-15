/*
 * Burgers_Plane_TS_ln_erk.hpp
 *
 *  Created on: 15 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#if SWEET_PARAREAL
#include <parareal/Parareal_GenericData.hpp>
#endif

class Burgers_Plane_TS_interface
{
public:
	virtual void run_timestep(
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables
			///PlaneData_Spectral &io_u_prev,	///< prognostic variables
			///PlaneData_Spectral &io_v_prev,	///< prognostic variables

			double i_fixed_dt,
			double i_simulation_timestamp
	) = 0;

#if SWEET_PARAREAL
	void run_timestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		PlaneData_Spectral u = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
		PlaneData_Spectral v = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]);

		run_timestep(u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = u;
		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = v;

	}

	// for parareal SL
	virtual void set_previous_solution(
				PlaneData_Spectral &i_u_prev,
				PlaneData_Spectral &i_v_prev
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			Parareal_GenericData* i_data
	)
	{
		PlaneData_Spectral u_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0];
		PlaneData_Spectral v_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1];

		set_previous_solution(u_prev, v_prev);
	};
#endif


};



#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_LN_ERK_HPP_ */
