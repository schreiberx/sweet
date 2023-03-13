/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"



class PDESWESphereTS_l_erk_n_erk	: public PDESWESphereTS_BaseInterface
{
public:
	sweet::TimesteppingExplicitRKSphereData timestepping_rk_linear;
	sweet::TimesteppingExplicitRKSphereData timestepping_rk_nonlinear;

	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup_main(
			sweet::SphereOperators *io_ops,
			int i_order,	///< order of RK time stepping method
			int i_order2
	);

	PDESWESphereTS_l_erk_n_erk();

	virtual ~PDESWESphereTS_l_erk_n_erk();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi_pert,
			sweet::SphereData_Spectral &io_vort,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

public:
	void euler_timestep_update_linear(
			const sweet::SphereData_Spectral &i_h,
			const sweet::SphereData_Spectral &i_u,
			const sweet::SphereData_Spectral &i_v,

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_nonlinear(
			const sweet::SphereData_Spectral &i_h,
			const sweet::SphereData_Spectral &i_u,
			const sweet::SphereData_Spectral &i_v,

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_nonlinear(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vort,
			sweet::SphereData_Spectral &io_div,

			double i_dt,
			double i_simulation_timestamp
	);

public:
	void euler_timestep_update(
			const sweet::SphereData_Spectral &i_phi,
			const sweet::SphereData_Spectral &i_vort,
			const sweet::SphereData_Spectral &i_div,

			sweet::SphereData_Spectral &o_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_vort_t,	///< time updates
			sweet::SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

};

#endif
