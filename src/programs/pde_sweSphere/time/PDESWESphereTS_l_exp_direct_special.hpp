#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_EXP_SPECIAL_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_EXP_SPECIAL_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_lg_exp_direct.hpp"
#include "PDESWESphereTS_ln_erk_split_vd.hpp"



class PDESWESphereTS_l_exp_direct_special	: public PDESWESphereTS_BaseInterface
{
	/*
	 * This class acts as a wrapper around the lg_exp_direct
	 * method and an ETDnRK version (lg + lc) to automatically
	 * choose the right one.
	 *
	 * Note, that the lc ETD method is actually not a direct method.
	 */

public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup_main(
			const sweet::SphereOperators *io_ops,
			int i_order,	///< order of RK time stepping method
			bool i_use_coriolis,		///< Include Coriolis term
			const std::string i_function_name	///< phi/ups function
	);

public:
	bool implementsTimesteppingMethod(
		const std::string &i_timestepping_method
	) override;

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override;

public:

	std::string getIDString() override
	{
		return "l_exp_special";
	}

	bool use_coriolis;

	PDESWESphereTS_lg_exp_direct timestepping_lg_exp_phi0;
	PDESWESphereTS_lg_exp_direct timestepping_lg_exp_phi1;
	PDESWESphereTS_lg_exp_direct timestepping_lg_exp_phi2;

	PDESWESphereTS_lg_exp_direct timestepping_lg_exp_ups1;
	PDESWESphereTS_lg_exp_direct timestepping_lg_exp_ups2;
	PDESWESphereTS_lg_exp_direct timestepping_lg_exp_ups3;

	PDESWESphereTS_ln_erk_split_vd timestepping_lc_erk;

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	void euler_timestep_store_update_lc(
			const sweet::SphereData_Spectral &i_phi_pert,
			const sweet::SphereData_Spectral &i_vrt,
			const sweet::SphereData_Spectral &i_div,
			sweet::SphereData_Spectral &o_phi_pert,
			sweet::SphereData_Spectral &o_vrt,
			sweet::SphereData_Spectral &o_div,
			double i_simulation_timestamp
	);

	PDESWESphereTS_l_exp_direct_special();

	virtual ~PDESWESphereTS_l_exp_direct_special();
};

#endif

