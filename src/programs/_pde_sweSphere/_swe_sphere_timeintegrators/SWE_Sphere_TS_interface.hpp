/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_INTERFACE_NEW_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_INTERFACE_NEW_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#if SWEET_PARAREAL || SWEET_XBRAID
#include <sweet/parareal/Parareal_GenericData.hpp>
#endif


class SWE_Sphere_TS_interface
{
// needed for parareal (instead of using directily shackDict.disc.timestepping_method)
protected:
	std::string timestepping_method;
	int timestepping_order;
	int timestepping_order2;

public:

	/*
	 * Automatic setup based on shackDict and operator
	 */
	virtual void setup_auto() = 0;

	/*
	 * Timestepping interface used by main timestepping loop
	 */
	virtual void run_timestep(
			sweet::SphereData_Spectral &io_h,	///< prognostic variables
			sweet::SphereData_Spectral &io_u,	///< prognostic variables
			sweet::SphereData_Spectral &io_v,	///< prognostic variables

			double i_fixed_dt,
			double i_simulation_timestamp
	) = 0;

	virtual bool implements_timestepping_method(
			const std::string &i_timestepping_method
		) = 0;

	virtual std::string string_id() = 0;

	virtual void print_help()
	{
	}

#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE)
	void run_timestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		sweet::SphereData_Spectral h = *(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[0]);
		sweet::SphereData_Spectral u = *(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[1]);
		sweet::SphereData_Spectral v = *(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[2]);

		run_timestep(h, u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[0]) = h;
		*(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[1]) = u;
		*(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[2]) = v;

	}



	// for parareal SL
	virtual void set_previous_solution(
				sweet::SphereData_Spectral &i_phi_prev,
				sweet::SphereData_Spectral &i_vrt_prev,
				sweet::SphereData_Spectral &i_div_prev
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			Parareal_GenericData* i_data
	)
	{
		sweet::SphereData_Spectral phi_prev = *i_data->get_pointer_to_data_SphereData_Spectral()->simfields[0];
		sweet::SphereData_Spectral vrt_prev = *i_data->get_pointer_to_data_SphereData_Spectral()->simfields[1];
		sweet::SphereData_Spectral div_prev = *i_data->get_pointer_to_data_SphereData_Spectral()->simfields[2];

		set_previous_solution(phi_prev, vrt_prev, div_prev);
	};
#endif



	virtual ~SWE_Sphere_TS_interface()
	{
	}

};

#endif
