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
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackDict.disc.timestepping_order;
		//timestepping_order2 = shackDict.disc.timestepping_order2;
		return i_timestepping_method == "lg_exp_lc_exp";
	}

public:
	std::string string_id()
	{
		return "lg_exp_lc_exp";
	}

	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &op;

	int timestepping_order;

	PDESWESphereTS_lg_exp_direct timestepping_lg_exp;
	PDESWESphereTS_ln_erk_split_vd timestepping_lc_erk;



public:
	PDESWESphereTS_lg_exp_lc_taylor(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void setup_auto();

	void run_timestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void run_timestep_lc(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void run_timestep_lg(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~PDESWESphereTS_lg_exp_lc_taylor();
};

#endif
