/*
 * PDESWESphereTS_l_erk_na_erk_vd.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_NA_ERK_VD_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_NA_ERK_VD_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_ln_erk_split_vd.hpp"



class PDESWESphereTS_l_erk_na_erk_vd	: public PDESWESphereTS_BaseInterface
{
public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	PDESWESphereTS_ln_erk_split_vd l_erk_split_vd;
	PDESWESphereTS_ln_erk_split_vd na_erk_split_vd;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		l_erk_split_vd.shackRegistration(io_shackDict);
		na_erk_split_vd.shackRegistration(io_shackDict);
		return true;
	}

public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup_main(
			sweet::SphereOperators *io_ops,
			int i_order,	///< order of RK time stepping method
			int i_order2
	);

public:
	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	PDESWESphereTS_l_erk_na_erk_vd();


	virtual ~PDESWESphereTS_l_erk_na_erk_vd();
};

#endif
