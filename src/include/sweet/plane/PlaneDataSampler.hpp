/*
 * Sampler2D.hpp
 *
 *  Created on: 4 Dec 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_
#define SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_

//#include "PlaneDataComplex.hpp"

class PlaneDataSampler
{
public:
	double domain_size[2];
	int res[2];
	PlaneDataConfig *planeDataConfig;

private:
	double scale_factor[2];


public:
	PlaneDataSampler(
		double i_domain_size[2],
		PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(i_planeDataConfig != nullptr);

		planeDataConfig = i_planeDataConfig;
		setup(i_domain_size, planeDataConfig);
	}


	PlaneDataSampler()
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
		PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(i_planeDataConfig != nullptr);
		planeDataConfig = i_planeDataConfig;

		domain_size[0] = i_domain_size[0];
		domain_size[1] = i_domain_size[1];

		res[0] = i_planeDataConfig->physical_res[0];
		res[1] = i_planeDataConfig->physical_res[1];

		scale_factor[0] = (double)i_planeDataConfig->physical_res[0] / i_domain_size[0];
		scale_factor[1] = (double)i_planeDataConfig->physical_res[1] / i_domain_size[1];
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
			PlaneData &i_u,				///< u data in C-grid
			PlaneData &i_v,				///< v data in C-grid
			PlaneData &o_u,				///< u data in A-grid
			PlaneData &o_v,				///< v data in A-grid
			int i_method=0
	)
	{
		// position of A grid points (h)
		PlaneData pos_x(i_u.planeDataConfig);
		PlaneData pos_y(i_u.planeDataConfig);

		//Initialise output
		o_u=i_u;
		o_v=i_v;

		assert(res[0] > 0);
		assert(scale_factor[0] > 0);


		pos_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = ((double)i+0.5)/scale_factor[0];
			});
		pos_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = ((double)j+0.5)/scale_factor[1];
			});


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
			PlaneData &i_data,				///< sampling data
			PlaneData &i_pos_x,				///< x positions of interpolation points
			PlaneData &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	            ///< sampling position
			PlaneData &o_data,				///< output values
			double i_shift_x = 0.0,            ///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(scale_factor[0] > 0);

		const std::size_t size = i_pos_x.planeDataConfig->physical_array_data_number_of_elements;

		i_data.request_data_physical();

		// iterate over all positions
//#pragma omp parallel for OPENMP_PAR_SIMD
#pragma omp parallel for
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
			double pos_x = wrapPeriodic(i_pos_x.physical_space_data[pos_idx]*scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.physical_space_data[pos_idx]*scale_factor[1] + i_shift_y, (double)res[1]);
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

				p[0] = i_data.physical_get(idx_j, idx_i[0]);
				p[1] = i_data.physical_get(idx_j, idx_i[1]);
				p[2] = i_data.physical_get(idx_j, idx_i[2]);
				p[3] = i_data.physical_get(idx_j, idx_i[3]);

				q[kj] = p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
				//std::cout << idx_j << " "<< q[kj] << std::endl;
				idx_j = wrapPeriodic(idx_j+1, res[1]);

			}
			//std::cout << q[0] << " " << q[1] <<  " " << q[2] <<  " " <<  q[3] <<  " " << y << std::endl;
			double value = q[1] + 0.5 * y*(q[2] - q[0] + y*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + y*(3.0*(q[1] - q[2]) + q[3] - q[0])));
			//std::cout << value << std::endl;
			//std::cout  << std::endl;
			o_data.physical_space_data[pos_idx] = value;
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		o_data.physical_space_data_valid = true;
		o_data.spectral_space_data_valid = false;
#endif
	}


public:
	void bilinear_scalar(
			PlaneData &i_data,				///< sampling data
			PlaneData &i_pos_x,				///< x positions of interpolation points
			PlaneData &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	///< sampling position
			PlaneData &o_data,				///< output values
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		/*
		 * SHIFT - important
		 * for C grid, to interpolate given u data, use i_shift_x = 0.0,  i_shift_y = -0.5
		 *             to interpolate given v data, use i_shift_y = -0.5, i_shift_y = 0.0
		 *  pay attention to the negative shift, which is necessary because the staggered grids are positively shifted
		 *  and this shift has to be removed for the interpolation
		 */


		i_data.request_data_physical();
		std::size_t size = i_pos_x.planeDataConfig->physical_array_data_number_of_elements;

		assert(size != 0);

		// iterate over all positions
//#pragma omp parallel for OPENMP_PAR_SIMD
#pragma omp parallel for
		for (std::size_t pos_idx = 0; pos_idx < size; pos_idx++)
		{
			// load position to interpolate
			double pos_x = wrapPeriodic(i_pos_x.physical_space_data[pos_idx]*scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.physical_space_data[pos_idx]*scale_factor[1] + i_shift_y, (double)res[0]);

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
				p[0] = i_data.physical_get(idx_j, idx_i[0]);
				p[1] = i_data.physical_get(idx_j, idx_i[1]);

				q[kj] = p[0] + x*(p[1]-p[0]);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
				//std::cout<< p[0] << p[1] << x << std::endl;
			}

			double value = q[0] + y*(q[1]-q[0]);

			o_data.physical_space_data[pos_idx] = value;
		}

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		o_data.physical_space_data_valid = true;
		o_data.spectral_space_data_valid = false;
#endif
	}


public:
	const PlaneData bilinear_scalar(
			PlaneData &i_data,				///< sampling data
			PlaneData &i_pos_x,				///< x positions of interpolation points
			PlaneData &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		PlaneData out(i_data.planeDataConfig);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}

public:
	const PlaneData bicubic_scalar(
			PlaneData &i_data,				///< sampling data
			PlaneData &i_pos_x,				///< x positions of interpolation points
			PlaneData &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		PlaneData out(planeDataConfig);
		bicubic_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}

#if 0
	/*
	 *
	 *  Bicubic interpolation routine
	 *  Receives complex arrays and returns complex cartesian data
	 *    with null imag part.
	 *   Data MUST be in cartesian space!!!!
	 */
public:
	const PlaneDataComplex bicubic_scalar(
			PlaneDataComplex &i_data,				///< sampling data
			PlaneData &i_pos_x,				///< x positions of interpolation points
			PlaneData &i_pos_y,				///< y positions of interpolation points
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{

		PlaneData data(planeDataConfig); //i_data in real data_array structure
		PlaneData out(planeDataConfig); // interpolated data in real data_array structure
		PlaneDataComplex out_cmp(planeDataConfig); // complex output of interpolated data

		// The data needs to be Cartesian space!!
		//data_cmp=i_data.toCart(); // do not use, since not secure

		//Put data into a Real PlaneData - called data
		i_data.toPlaneDatas_Real(data);

		// Do the interpolation
		bicubic_scalar(data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);

		//Convert back to complex array
		out_cmp.loadRealFromPlaneData(out);

		return out_cmp;
	}
#endif
};




#endif /* SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_ */
