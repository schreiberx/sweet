/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_LC_N_ETD_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_LC_N_ETD_UV_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_exp.hpp"
#include "PDESWESphereTS_ln_erk_split_uv.hpp"


class PDESWESphereTS_lg_exp_lc_n_etd_uv	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(sweet::SphereOperators *io_ops);

	bool setup(
			sweet::SphereOperators *io_ops,
			sweet::ShackExpIntegration *i_shackExpIntegration,
			int i_timestepping_order,
			double i_timestep_size,
			bool i_with_na,
			bool i_with_nr
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method);
	std::string getIDString();
	void printHelp();


	PDESWESphereTS_ln_erk_split_uv ts_ln_erk_split_uv;

	sweet::SphereData_Spectral NU_phi_prev, NU_vrt_prev, NU_div_prev;
	sweet::SphereData_Spectral NU_phi_prev_2, NU_vrt_prev_2, NU_div_prev_2;

	bool with_na;
	bool with_nr;

	PDESWESphereTS_l_exp ts_phi0_exp;
	PDESWESphereTS_l_exp ts_phi1_exp;
	PDESWESphereTS_l_exp ts_phi2_exp;
	PDESWESphereTS_l_exp ts_phi3_exp;

private:
	void euler_timestep_update_linear(
			const sweet::SphereData_Spectral &i_h,	///< prognostic variables
			const sweet::SphereData_Spectral &i_u,	///< prognostic variables
			const sweet::SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);



private:
	void euler_timestep_update_nonlinear(
			const sweet::SphereData_Spectral &i_h,	///< prognostic variables
			const sweet::SphereData_Spectral &i_u,	///< prognostic variables
			const sweet::SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);


public:
	PDESWESphereTS_lg_exp_lc_n_etd_uv();


	void runTimestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~PDESWESphereTS_lg_exp_lc_n_etd_uv();
};

#endif
