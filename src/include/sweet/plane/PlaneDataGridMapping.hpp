/*
 * PlaneDataGridMapping.hpp
 *
 *  Created on: 18 Jul 2017
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_PLANEDATAGRIDMAPPING_HPP_
#define SRC_INCLUDE_SWEET_PLANE_PLANEDATAGRIDMAPPING_HPP_

#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/ScalarDataArray.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneStaggering.hpp>

class PlaneDataGridMapping
{
public:
	//(x,y) grid points, refers to lower left corner of cells
	ScalarDataArray pos_ll_x, pos_ll_y;

	// Interpolation stuff
	PlaneDataSampler sampler2D;

	// Staggering
	Staggering staggering;

public:
	PlaneDataGridMapping()
	{
	}


	void setup(
			SimulationVariables i_simVars,
			PlaneDataConfig *i_planeDataConfig
	)
	{
		// ll  refers to lower left corner of the cell.
		pos_ll_x.setup(i_planeDataConfig->physical_array_data_number_of_elements);
		pos_ll_y.setup(i_planeDataConfig->physical_array_data_number_of_elements);

		std::size_t idx = 0;
		for (std::size_t j = 0; j < i_planeDataConfig->physical_res[1]; j++)
		{
			for (std::size_t i = 0; i < i_planeDataConfig->physical_res[0]; i++)
			{
				pos_ll_x.scalar_data[idx] = ((double)i)*i_simVars.sim.plane_domain_size[0]/(double)i_simVars.disc.space_res_physical[0];
				pos_ll_y.scalar_data[idx] = ((double)j)*i_simVars.sim.plane_domain_size[1]/(double)i_simVars.disc.space_res_physical[1];
				idx++;
			}
		}

		// Setup sampler for future interpolations
		sampler2D.setup(i_simVars.sim.plane_domain_size, i_planeDataConfig);

		if (i_simVars.disc.space_grid_use_c_staggering)
			staggering.setup_c_staggering();
		else
			staggering.setup_a_staggering();
	}


	void mapCtoA_u(
			const PlaneData_Spectral &i_src,
			PlaneData_Spectral &o_dst
	)
	{
		// remap solution to A grid
		sampler2D.bicubic_scalar(i_src.toPhys(), pos_ll_x, pos_ll_y, o_dst, staggering.u[0], staggering.u[1]);
	}


	void mapCtoA_v(
			const PlaneData_Spectral &i_src,
			PlaneData_Spectral &o_dst
	)
	{
		// remap solution to A grid
		sampler2D.bicubic_scalar(i_src.toPhys(), pos_ll_x, pos_ll_y, o_dst, staggering.v[0], staggering.v[1]);
	}




	void mapAtoC_u(
			const PlaneData_Spectral &i_src,
			PlaneData_Spectral &o_dst
	)
	{
		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src.toPhys(), pos_ll_x, pos_ll_y, o_dst, -staggering.u[0], -staggering.u[1]);
	}


	void mapAtoC_v(
			const PlaneData_Spectral &i_src,
			PlaneData_Spectral &o_dst
	)
	{
		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src.toPhys(), pos_ll_x, pos_ll_y, o_dst, -staggering.v[0], -staggering.v[1]);
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANE_PLANEDATAGRIDMAPPING_HPP_ */
