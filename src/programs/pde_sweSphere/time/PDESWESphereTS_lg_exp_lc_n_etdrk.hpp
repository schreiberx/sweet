/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_PDESWESphereTS_lg_rexi_lf_n_etdrk_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_PDESWESphereTS_lg_rexi_lf_n_etdrk_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_exp.hpp"
#include "PDESWESphereTS_lg_erk_lc_n_erk.hpp"


class PDESWESphereTS_lg_exp_lc_n_etdrk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup_main(
		const sweet::SphereOperators *io_ops,
		sweet::ShackExpIntegration *i_shackExpIntegration,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestepSize
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;


private:
	PDESWESphereTS_lg_erk_lc_n_erk ts_lg_erk_lc_n_erk;

	PDESWESphereTS_l_exp ts_phi0_rexi;
	PDESWESphereTS_l_exp ts_phi1_rexi;
	PDESWESphereTS_l_exp ts_phi2_rexi;

	PDESWESphereTS_l_exp ts_ups0_rexi;
	PDESWESphereTS_l_exp ts_ups1_rexi;
	PDESWESphereTS_l_exp ts_ups2_rexi;
	PDESWESphereTS_l_exp ts_ups3_rexi;


public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		ts_lg_erk_lc_n_erk.shackRegistration(io_shackDict);

		ts_phi0_rexi.shackRegistration(io_shackDict);
		ts_phi1_rexi.shackRegistration(io_shackDict);
		ts_phi2_rexi.shackRegistration(io_shackDict);

		ts_ups0_rexi.shackRegistration(io_shackDict);
		ts_ups1_rexi.shackRegistration(io_shackDict);
		ts_ups2_rexi.shackRegistration(io_shackDict);
		ts_ups3_rexi.shackRegistration(io_shackDict);
		return true;
	}


private:
	void euler_timestep_update_linear(
			const sweet::SphereData_Spectral &i_h,
			const sweet::SphereData_Spectral &i_u,
			const sweet::SphereData_Spectral &i_v,

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);



private:
	void euler_timestep_update_nonlinear(
			const sweet::SphereData_Spectral &i_h,
			const sweet::SphereData_Spectral &i_u,
			const sweet::SphereData_Spectral &i_v,

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_max_timestamp
	);


public:
	PDESWESphereTS_lg_exp_lc_n_etdrk();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	virtual ~PDESWESphereTS_lg_exp_lc_n_etdrk();
};

#endif
