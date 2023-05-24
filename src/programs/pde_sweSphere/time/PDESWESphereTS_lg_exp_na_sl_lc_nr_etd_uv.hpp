/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_NA_SL_LC_N_ETD_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_NA_SL_LC_N_ETD_UV_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include <sweet/core/time/ShackTimesteppingSemiLagrangianSphereData.hpp>
#include <sweet/core/time/TimesteppingSemiLagrangianSphereData.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_exp.hpp"
#include "PDESWESphereTS_ln_erk_split_uv.hpp"


class PDESWESphereTS_lg_exp_na_sl_lc_nr_etd_uv	: public PDESWESphereTS_BaseInterface
{
public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;
	void printHelp() override;

private:
	PDESWESphereTS_ln_erk_split_uv ts_ln_erk_split_uv;

	std::string timestepping_method = "";

private:
	enum NLRemainderTreatment_enum{
		NL_REMAINDER_IGNORE,
		NL_REMAINDER_NONLINEAR,
	};

	NLRemainderTreatment_enum nonlinear_remainder_treatment;

public:
	sweet::SphereData_Spectral U_phi_prev, U_vrt_prev, U_div_prev;

	PDESWESphereTS_l_exp ts_phi0_exp;
	PDESWESphereTS_l_exp ts_phi1_exp;
	PDESWESphereTS_l_exp ts_phi2_exp;


public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		ts_phi0_exp.shackRegistration(io_shackDict);
		ts_phi1_exp.shackRegistration(io_shackDict);
		ts_phi2_exp.shackRegistration(io_shackDict);
		return true;
	}



	sweet::TimesteppingSemiLagrangianSphereData semiLagrangian;

public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup(
			sweet::SphereOperators *io_ops,
			sweet::ShackExpIntegration *i_shackExpIntegration,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size,

			NLRemainderTreatment_enum i_nonlinear_remainder_treatment
	);


public:
	PDESWESphereTS_lg_exp_na_sl_lc_nr_etd_uv();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;

#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE)
	void set_previous_solution(
				sweet::SphereData_Spectral &i_phi_prev,
				sweet::SphereData_Spectral &i_vrt_prev,
				sweet::SphereData_Spectral &i_div_prev
	) override
	{
		U_phi_prev = i_phi_prev;
		U_vrt_prev = i_vrt_prev;
		U_div_prev = i_div_prev;
	}
#endif

	virtual ~PDESWESphereTS_lg_exp_na_sl_lc_nr_etd_uv();
};

#endif
