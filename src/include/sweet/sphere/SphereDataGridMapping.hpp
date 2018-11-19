/*
 * SphereDataGridMapping.hpp
 *
 *  Created on: 18 Jul 2017
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_SPHEREDATAGRIDMAPPING_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_SPHEREDATAGRIDMAPPING_HPP_

#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/ScalarDataArray.hpp>
#include <sweet/sphere/SphereStaggering.hpp>
#include <sweet/sphere/SphereDataSampler.hpp>

class SphereDataGridMapping
{
public:
	//(x,y) grid points, refers to lower left corner of cells
	ScalarDataArray pos_ll_x, pos_ll_y;

	// Interpolation stuff
	SphereDataSampler sampler2D;

	// Staggering
	SphereStaggering staggering;


public:
	SphereDataGridMapping()
	{
	}


	void setup(
			SimulationVariables i_simVars,
			SphereDataConfig *i_sphereDataConfig
	)
	{
		// ll  refers to lower left corner of the cell.
		pos_ll_x.setup(i_sphereDataConfig->physical_array_data_number_of_elements);
		pos_ll_y.setup(i_sphereDataConfig->physical_array_data_number_of_elements);

		std::size_t idx = 0;
		for (std::size_t j = 0; j < i_sphereDataConfig->physical_res[1]; j++)
		{
			for (std::size_t i = 0; i < i_sphereDataConfig->physical_res[0]; i++)
			{
				pos_ll_x.scalar_data[idx] = ((double)i)*i_simVars.sim.plane_domain_size[0]/(double)i_simVars.disc.space_res_physical[0];
				pos_ll_y.scalar_data[idx] = ((double)j)*i_simVars.sim.plane_domain_size[1]/(double)i_simVars.disc.space_res_physical[1];
				idx++;
			}
		}

		// Setup sampler for future interpolations
		sampler2D.setup(i_simVars.sim.plane_domain_size, i_sphereDataConfig);

		if (i_simVars.disc.space_grid_use_c_staggering)
			staggering.setup_c_staggering();
		else
			staggering.setup_a_staggering();
	}


	void mapCtoA_u(
			const SphereData &i_src,
			SphereData &o_dst
	)
	{
		// remap solution to A grid
		sampler2D.bicubic_scalar(i_src, pos_ll_x, pos_ll_y, o_dst, staggering.u[0], staggering.u[1]);
	}


	void mapCtoA_v(
			const SphereData &i_src,
			SphereData &o_dst
	)
	{
		// remap solution to A grid
		sampler2D.bicubic_scalar(i_src, pos_ll_x, pos_ll_y, o_dst, staggering.v[0], staggering.v[1]);
	}




	void mapAtoC_u(
			const SphereData &i_src,
			SphereData &o_dst
	)
	{
		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src, pos_ll_x, pos_ll_y, o_dst, -staggering.u[0], -staggering.u[1]);
	}


	void mapAtoC_v(
			const SphereData &i_src,
			SphereData &o_dst
	)
	{
		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src, pos_ll_x, pos_ll_y, o_dst, -staggering.v[0], -staggering.v[1]);
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_SPHEREDATAGRIDMAPPING_HPP_ */
