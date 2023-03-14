#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_EXP_LC_EXP_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_EXP_LC_EXP_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_lg_exp_direct.hpp"
#include "PDESWESphereTS_ln_erk_split_vd.hpp"



class PDESWESphereTS_lg_exp_lc_taylor	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup_main(
			sweet::SphereOperators *io_ops,
			int i_order	///< order of RK time stepping method
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		//timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		return i_timestepping_method == "lg_exp_lc_exp";
	}

public:
	std::string getIDString() override
	{
		return "lg_exp_lc_exp";
	}

	PDESWESphereTS_lg_exp_direct timestepping_lg_exp;
	PDESWESphereTS_ln_erk_split_vd timestepping_lc_erk;


public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_lg_exp.shackRegistration(io_shackDict);
		timestepping_lc_erk.shackRegistration(io_shackDict);
		return true;
	}


public:
	PDESWESphereTS_lg_exp_lc_taylor();

	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	void run_timestep_lc(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void run_timestep_lg(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~PDESWESphereTS_lg_exp_lc_taylor();
};

#endif
