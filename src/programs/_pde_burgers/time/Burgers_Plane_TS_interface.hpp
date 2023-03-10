/*
 * Burgers_Plane_TS_ln_erk.hpp
 *
 *  Created on: 15 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_INTERFACE_HPP_

#include <limits>
#include <sweet/core/plane/sweet::PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#if SWEET_PARAREAL || SWEET_XBRAID
#include <sweet/parareal/Parareal_GenericData.hpp>
#endif

class Burgers_Plane_TS_interface
{
public:
	virtual void runTimestep(
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables
			///sweet::PlaneData_Spectral &io_u_prev,	///< prognostic variables
			///sweet::PlaneData_Spectral &io_v_prev,	///< prognostic variables

			double i_fixed_dt,
			double i_simulation_timestamp
	) = 0;

#if (SWEET_PARAREAL && SWEET_PARAREAL_PLANE) || (SWEET_XBRAID && SWEET_XBRAID_PLANE)
	void runTimestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		sweet::PlaneData_Spectral u = *(io_data->get_pointer_to_data_sweet::PlaneData_Spectral()->simfields[0]);
		sweet::PlaneData_Spectral v = *(io_data->get_pointer_to_data_sweet::PlaneData_Spectral()->simfields[1]);

		runTimestep(u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_sweet::PlaneData_Spectral()->simfields[0]) = u;
		*(io_data->get_pointer_to_data_sweet::PlaneData_Spectral()->simfields[1]) = v;

	}

	// for parareal SL
	virtual void set_previous_solution(
				sweet::PlaneData_Spectral &i_u_prev,
				sweet::PlaneData_Spectral &i_v_prev
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			Parareal_GenericData* i_data
	)
	{
		sweet::PlaneData_Spectral u_prev = *i_data->get_pointer_to_data_sweet::PlaneData_Spectral()->simfields[0];
		sweet::PlaneData_Spectral v_prev = *i_data->get_pointer_to_data_sweet::PlaneData_Spectral()->simfields[1];

		set_previous_solution(u_prev, v_prev);
	};
#endif


};



#endif
