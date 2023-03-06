/*
 * PDESWESphereTS_l_erk_na_erk_vd.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
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
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();


	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &ops;

	int timestepping_order;
	int timestepping_order2;

	PDESWESphereTS_ln_erk_split_vd *l_erk_split_vd = nullptr;
	PDESWESphereTS_ln_erk_split_vd *na_erk_split_vd = nullptr;

public:
	PDESWESphereTS_l_erk_na_erk_vd(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_order2
	);


	void setup_auto();

	void run_timestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~PDESWESphereTS_l_erk_na_erk_vd();
};

#endif
