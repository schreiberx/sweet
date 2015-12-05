/*
 * Sampler2D.hpp
 *
 *  Created on: 4 Dec 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_SAMPLER2D_HPP_
#define SRC_INCLUDE_SWEET_SAMPLER2D_HPP_


class Sampler2D
{
	double domain_size[2];
	int res[2];

	double scale_factor[2];


public:
	Sampler2D(
		double i_domain_size[2],
		std::size_t i_res[2]
	)
	{
		setup(i_domain_size, i_res);
	}

	Sampler2D()
	{
	}


public:
	void setup(
		double i_domain_size[2],
		std::size_t i_res[2]
	)
	{
		domain_size[0] = i_domain_size[0];
		domain_size[1] = i_domain_size[1];

		res[0] = i_res[0];
		res[1] = i_res[1];

		scale_factor[0] = (double)i_res[0] / i_domain_size[0];
		scale_factor[1] = (double)i_res[1] / i_domain_size[1];
	}

private:
	/**
	 * wrap the position i in a periodic domain of size d
	 */
	inline
	int wrapPeriodic(int i, int i_res)
	{
		while (i < 0)
			i += i_res;
		while (i >= i_res)
			i -= i_res;

		assert(i >= 0 && i < i_res);
		return i;
	}


public:
	void bicubic_scalar(
			DataArray<2> &i_data,				///< sampling data
			DataArray<2>* i_pos[2],	            ///< sampling position
			DataArray<2> &o_data,				///< output values
			double i_shift_x = 0.0,            ///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		// iterate over all positions
#pragma omp parallel for OPENMP_SIMD
		for (int pos_idx = 0; pos_idx < i_pos[0]->resolution[0]*i_pos[0]->resolution[1]; pos_idx++)
		{
			// load position to interpolate
			double pos_x = i_pos[0]->array_data_cartesian_space[pos_idx]*scale_factor[0] + i_shift_x;
			double pos_y = i_pos[1]->array_data_cartesian_space[pos_idx]*scale_factor[1] + i_shift_y;

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[4];
			{
				int i = (int)pos_x-1;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[2] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[3] = wrapPeriodic(i, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = wrapPeriodic((int)pos_y-1, res[1]);

			double q[4];
			for (int kj = 0; kj < 4; kj++)
			{
				double p[4];
				p[0] = i_data.get(idx_j, idx_i[0]);
				p[1] = i_data.get(idx_j, idx_i[1]);
				p[2] = i_data.get(idx_j, idx_i[2]);
				p[3] = i_data.get(idx_j, idx_i[3]);

				q[kj] = p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double value = q[1] + 0.5 * y*(q[2] - q[0] + y*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + y*(3.0*(q[1] - q[2]) + q[3] - q[0])));
			o_data.array_data_cartesian_space[pos_idx] = value;
		}
	}

public:
	void bilinear_scalar(
			DataArray<2> &i_data,				///< sampling data
			DataArray<2>* i_pos[2],	///< sampling position
			DataArray<2> &o_data,				///< output values
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		// iterate over all positions
#pragma omp parallel for OPENMP_SIMD
		for (int pos_idx = 0; pos_idx < i_pos[0]->resolution[0]*i_pos[0]->resolution[1]; pos_idx++)
		{
			// load position to interpolate
			double pos_x = i_pos[0]->array_data_cartesian_space[pos_idx]*scale_factor[0] + i_shift_x;
			double pos_y = i_pos[1]->array_data_cartesian_space[pos_idx]*scale_factor[1] + i_shift_y;

			/**
			 * For the interpolation, we assume node-aligned values
			 *
			 * x0  x1  x2  x3
			 * |---|---|---|---
			 * 0   2   4   6    <- positions and associated values e.g. for domain size 8
			 */

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 4 times
			int idx_i[2];
			{
				int i = (int)pos_x-1;

				i = wrapPeriodic(i, res[0]);
				idx_i[0] = wrapPeriodic(i, res[0]);

				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);

			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = wrapPeriodic((int)pos_y-1, res[1]);

			double q[2];
			for (int kj = 0; kj < 4; kj++)
			{
				double p[2];
				p[0] = i_data.get(idx_j, idx_i[0]);
				p[1] = i_data.get(idx_j, idx_i[1]);

				q[kj] = p[1] + x*p[2];

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}

			double value = q[1] + y*q[2];
			o_data.array_data_cartesian_space[pos_idx] = value;
		}
	}

};
#endif /* SRC_INCLUDE_SWEET_SAMPLER2D_HPP_ */
