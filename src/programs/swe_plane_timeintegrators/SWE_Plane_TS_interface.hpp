/*
 * SWE_Plane_TS_interface.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#if SWEET_PARAREAL
#include <parareal/Parareal_GenericData.hpp>
#endif

class SWE_Plane_TS_interface
{
public:
	virtual void run_timestep(
			PlaneData_Spectral &io_h_pert,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,		///< time step size
			double i_sim_timestamp
	) = 0;

#if SWEET_PARAREAL
	void run_timestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		PlaneData_Spectral h_pert = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
		PlaneData_Spectral u = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]);
		PlaneData_Spectral v = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2]);

		run_timestep(h_pert, u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = h_pert;
		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = u;
		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2]) = v;

	}

	// for parareal SL
	virtual void set_previous_solution(
				PlaneData_Spectral &i_h_prev,
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
		PlaneData_Spectral h_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0];
		PlaneData_Spectral u_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1];
		PlaneData_Spectral v_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2];

		set_previous_solution(h_prev, u_prev, v_prev);
	};
#endif

};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_INTERFACE_HPP_ */
