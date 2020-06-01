/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_INTERFACE_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_INTERFACE_HPP_

#include <sweet/ScalarDataArray.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>


class BenchmarksSphereAdvection_interface
{
public:
	virtual void setup(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
		) = 0;

	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;

	virtual void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{
		SWEETError("Not implemented for this benchmark");
	}

	virtual void get_initial_state(
		std::vector<SphereData_Spectral*> &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{
		SWEETError("Not implemented for this benchmark");
	}

	virtual std::string get_help() = 0;


	virtual void get_reference_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div,
		double i_timestamp
	)
	{
		SWEETError("Not implemented for this benchmark");
	}

	/*
	 * Return number of prognostic fields to be used
	 */
	virtual int get_num_prognostic_fields()
	{
		return 1;
	}

	virtual bool has_time_varying_state()
	{
		return false;
	}

	virtual void sl_compute_departure_3rd_order(
			const ScalarDataArray &i_pos_lon_A,	///< longitude coordinate to compute the velocity for
			const ScalarDataArray &i_pos_lat_A,	///< latitude coordinate to compute the velocity for
			ScalarDataArray &o_pos_lon_D,		///< velocity along longitude
			ScalarDataArray &o_pos_lat_D,		///< velocity along latitude
			double i_dt,
			double i_timestamp_arrival			///< timestamp at arrival point
	)
	{
		SWEETError("Not implemented for this benchmark");
	}


	virtual ~BenchmarksSphereAdvection_interface()
	{
	}
};

#endif
