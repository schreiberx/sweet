/*
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_NA_ERK_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_NA_ERK_UV_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_ln_erk_split_uv.hpp"



class PDESWESphereTS_l_erk_na_erk_uv	: public PDESWESphereTS_BaseInterface
{
public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method);
	std::string getIDString();

	PDESWESphereTS_ln_erk_split_uv *l_erk_split_uv = nullptr;
	PDESWESphereTS_ln_erk_split_uv *na_erk_split_uv = nullptr;

public:
	bool setup_auto(sweet::SphereOperators *io_ops);

	bool setup(
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
	);

	void clear();

	PDESWESphereTS_l_erk_na_erk_uv();

	virtual ~PDESWESphereTS_l_erk_na_erk_uv();
};

#endif
