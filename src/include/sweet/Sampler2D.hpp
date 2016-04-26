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
public:
	double domain_size[2];
	int res[2];

private:
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
		res[0] = -1;
		res[1] = -1;

		scale_factor[0] = -1;
		scale_factor[1] = -1;

		domain_size[0] = -1;
		domain_size[1] = -1;
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

public:
	/**
	 * wrap the position i in a periodic domain of size i_res
	 */
	template <typename T>
	inline
	T wrapPeriodic(T i, T i_res)
	{
		while (i < 0)
		{
			i += i_res;
//			std::cout << i << std::endl;
		}
		while (i >= i_res)
			i -= i_res;

		assert(i >= 0 && i < i_res);

		return i;
	}


public:
	// NOT NEEDED AND PRO, just use directly the interpolation rotines with shifts...
	void remap_gridC2A(
			DataArray<2> &i_u,				///< u data in C-grid
			DataArray<2> &i_v,				///< v data in C-grid
			DataArray<2> &o_u,				///< u data in A-grid
			DataArray<2> &o_v,				///< v data in A-grid
			int i_method=0
	)
	{
		// position of A grid points (h)
		DataArray<2> pos_x(i_u.resolution);
		DataArray<2> pos_y(i_u.resolution);

		//Initialise output
		o_u=i_u;
		o_v=i_v;

		assert(res[0] > 0);
		assert(scale_factor[0] > 0);

		/* std::cout<< "test " << res[0] << scale_factor[0] <<std::endl;
		* std::cout<< i_u <<std::endl;
		* std::cout<<std::endl;
		* std::cout<< i_v <<std::endl;
		* std::cout<<std::endl;
		*/
		for (int j = 0; j < res[1]; j++)
		{
			for (int i = 0; i < res[0]; i++)
			{
		    	// h position - A grid
				pos_x.set(j, i, ((double)i+0.5)/scale_factor[0]); //*simVars.sim.domain_size[0];
				pos_y.set(j, i, ((double)j+0.5)/scale_factor[1]); //*simVars.sim.domain_size[1];
				//std::cout<< "i " << i << " j " << j << " x " << x <<" y "<< y <<std::endl;
			}
		}


		if(i_method>0){
			//Do bicubic interpolation
			bicubic_scalar(i_u, pos_x, pos_y, o_u, 0.0, -0.5);
			bicubic_scalar(i_v, pos_x, pos_y, o_v, -0.5, 0.0);
		}
		else{
			bilinear_scalar(i_u, pos_x, pos_y, o_u, 0.0, -0.5);
			bilinear_scalar(i_v, pos_x, pos_y, o_v, -0.5, 0.0);
		}

		//std::cout<< o_u <<std::endl;
		//std::cout<< o_v <<std::endl;
	}


public:
	void bicubic_scalar(
			DataArray<2> &i_data,				///< sampling data
			DataArray<2> &i_pos_x,				///< x positions of interpolation points
			DataArray<2> &i_pos_y,				///< y positions of interpolation points
			//DataArray<2>* i_pos[2],	            ///< sampling position
			DataArray<2> &o_data,				///< output values
			double i_shift_x = 0.0,            ///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(scale_factor[0] > 0);

		const std::size_t size = i_pos_x.resolution[0]*i_pos_x.resolution[1];

		// iterate over all positions
#pragma omp parallel for OPENMP_SIMD
		for (std::size_t pos_idx = 0; pos_idx < size; pos_idx++)
		{
			/*
			 * load position to interpolate
			 * posx stores all x-coordinates of the arrival points
			 * posy stores all y-coordinates of the arrival points
			 *
			 * Both are arrays and matching array indices (pos_idx) below index the coordinates for the same point.
			 *
			 * Scale factor (Nx/dx, Ny/dy) maps from the physical space to the array space.
			 * The array space is from [0; N[
			 *
			 * shift_x/y is operating in array space. Hence, staggered grid can be
			 * implemented by setting this to 0.5 or -0.5
			 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
			 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
			 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
			 *  and this shift has to be removed for the interpolation
			 */
			//double pos_x = i_pos[0]->array_data_cartesian_space[pos_idx]*scale_factor[0] + i_shift_x;
			//double pos_y = i_pos[1]->array_data_cartesian_space[pos_idx]*scale_factor[1] + i_shift_y;
			double pos_x = wrapPeriodic(i_pos_x.array_data_cartesian_space[pos_idx]*scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.array_data_cartesian_space[pos_idx]*scale_factor[1] + i_shift_y, (double)res[1]);
			//std::cout << pos_idx << " x " << pos_x << " y " << pos_y << " res " << res[1] << std::endl;
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
			//std::cout << idx_i[0] << idx_i[1] <<idx_i[2] << idx_i[3]<< std::endl;
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
				//std::cout << idx_j << " "<< q[kj] << std::endl;
				idx_j = wrapPeriodic(idx_j+1, res[1]);

			}
			//std::cout << q[0] << " " << q[1] <<  " " << q[2] <<  " " <<  q[3] <<  " " << y << std::endl;
			double value = q[1] + 0.5 * y*(q[2] - q[0] + y*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + y*(3.0*(q[1] - q[2]) + q[3] - q[0])));
			//std::cout << value << std::endl;
			//std::cout  << std::endl;
			o_data.array_data_cartesian_space[pos_idx] = value;
		}

#if SWEET_USE_SPECTRAL_SPACE
		o_data.array_data_cartesian_space_valid = true;
		o_data.array_data_spectral_space_valid = false;
#endif
	}


public:
	void bilinear_scalar(
			DataArray<2> &i_data,				///< sampling data
			DataArray<2> &i_pos_x,				///< x positions of interpolation points
			DataArray<2> &i_pos_y,				///< y positions of interpolation points
			//DataArray<2>* i_pos[2],	///< sampling position
			DataArray<2> &o_data,				///< output values
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		/* SHIFT - important
		 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
		 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
		 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
		 *  and this shift has to be removed for the interpolation
		 */


		const std::size_t size = i_pos_x.resolution[0]*i_pos_x.resolution[1];

		assert(size != 0);

		// iterate over all positions
#pragma omp parallel for OPENMP_SIMD
		for (std::size_t pos_idx = 0; pos_idx < size; pos_idx++)
		{
			// load position to interpolate
			double pos_x = wrapPeriodic(i_pos_x.array_data_cartesian_space[pos_idx]*scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.array_data_cartesian_space[pos_idx]*scale_factor[1] + i_shift_y, (double)res[0]);

			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */
			// compute x/y position
			double x = pos_x - floor(pos_x);
			double y = pos_y - floor(pos_y);

			// precompute x-position indices since they are reused 2 times
			int idx_i[2];
			{
				int i = (int)pos_x;
				idx_i[0] = i;

				i = wrapPeriodic(i+1, res[0]);
				idx_i[1] = wrapPeriodic(i, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = pos_y;//wrapPeriodic((int)pos_y, res[1]);

			double q[2];
			for (int kj = 0; kj < 2; kj++)
			{
				double p[2];
				p[0] = i_data.get(idx_j, idx_i[0]);
				p[1] = i_data.get(idx_j, idx_i[1]);

				q[kj] = p[0] + x*(p[1]-p[0]);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
				//std::cout<< p[0] << p[1] << x << std::endl;
			}

			double value = q[0] + y*(q[1]-q[0]);

			o_data.array_data_cartesian_space[pos_idx] = value;
		}

#if SWEET_USE_SPECTRAL_SPACE
		o_data.array_data_cartesian_space_valid = true;
		o_data.array_data_spectral_space_valid = false;
#endif
	}


public:
	const DataArray<2> bilinear_scalar(
			DataArray<2> &i_data,				///< sampling data
			DataArray<2> &i_pos_x,				///< x positions of interpolation points
			DataArray<2> &i_pos_y,				///< y positions of interpolation points
			//DataArray<2>* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		DataArray<2> out(i_data.resolution);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}

};
#endif /* SRC_INCLUDE_SWEET_SAMPLER2D_HPP_ */
