/*
 * Sampler2D.hpp
 *
 *  Created on: 4 Dec 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_
#define SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_

#include <sweet/ScalarDataArray.hpp>
//#include "PlaneDataComplex.hpp"


/**
 * this is a sampler class which provides method to provide
 * interpolated sampled values on 2D physical data which
 * is provided by PlaneData
 */
class PlaneDataSampler
{
public:
	double domain_size[2];			/// real physical size of the domain
	int res[2];						/// resolution of domain
	const PlaneDataConfig *planeDataConfig;

private:
	double cached_scale_factor[2];			/// cached parameters for sampling


public:
	PlaneDataSampler(
		double i_domain_size[2],	/// real physical size of the domain
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(i_planeDataConfig != nullptr);

		planeDataConfig = i_planeDataConfig;
		setup(i_domain_size, planeDataConfig);
	}


	PlaneDataSampler()
	{
		planeDataConfig = nullptr;

		res[0] = -1;
		res[1] = -1;

		cached_scale_factor[0] = -1;
		cached_scale_factor[1] = -1;

		domain_size[0] = -1;
		domain_size[1] = -1;
	}



public:
	void setup(
		double i_domain_size[2],	/// real physical size of the domain
		const PlaneDataConfig *i_planeDataConfig
	)
	{
		assert(i_planeDataConfig != nullptr);
		planeDataConfig = i_planeDataConfig;

		domain_size[0] = i_domain_size[0];
		domain_size[1] = i_domain_size[1];

		res[0] = i_planeDataConfig->physical_res[0];
		res[1] = i_planeDataConfig->physical_res[1];

		cached_scale_factor[0] = (double)i_planeDataConfig->physical_res[0] / i_domain_size[0];
		cached_scale_factor[1] = (double)i_planeDataConfig->physical_res[1] / i_domain_size[1];
	}

public:
	/**
	 * wrap the position i in a periodic domain of size i_res
	 */
#if 0

#error "This modulo operation doesn't work!"
	inline
	int wrapPeriodic(int i, int i_res)
	{
		return (i + i_res*10) % i_res;
	}

	inline
	double wrapPeriodic(double i, double i_res)
	{
		return fmodf(i + i_res*10.0f, i_res);
	}

#else

	template <typename T>
	inline
	double wrapPeriodic(T i, T i_res)
	{
#if 1
		int c = 10;
		while (i < 0 && c-- > 0)
			i += i_res;

		int d = 10;
		while (i >= i_res && d-- > 0)
			i -= i_res;
#elif 1

			i = (i + i_res*10) % i_res;
		else if (typeid(T) == typeid(double))
			i = fmod(i + i_res*10.0, i_res);
		else
			i = fmodf(i + i_res*10.0f, i_res);

#else
		if (i < 0)
			i += i_res;

		if (i >= i_res)
			i -= i_res;
#endif

		if (i < 0 || i >= i_res)
			FatalError("Stopping here: Probably an unstable velocity field since more than one periodic movement exists.");

		assert(i >= 0 && i < i_res);

		return i;
	}
#endif


public:
	void bicubic_scalar(
			const PlaneData &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		i_data.request_data_physical();

		std::size_t max_pos_idx = i_pos_x.number_of_elements;

#if SWEET_DEBUG
#if SWEET_SPACE_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM

		if (omp_get_num_threads() > 1)
			FatalError("Are we already in parallel region? Threading race conditions likely!");
#endif
#endif

		// iterate over all positions in parallel
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t pos_idx = 0; pos_idx < max_pos_idx; pos_idx++)
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
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

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

// TODO: Each 2nd wrapPeriodic is obsolete!!!
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

				p[0] = i_data.p_physical_get(idx_j, idx_i[0]);
				p[1] = i_data.p_physical_get(idx_j, idx_i[1]);
				p[2] = i_data.p_physical_get(idx_j, idx_i[2]);
				p[3] = i_data.p_physical_get(idx_j, idx_i[3]);

				q[kj] = p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

				idx_j = wrapPeriodic(idx_j+1, res[1]);
			}
			double value = q[1] + 0.5 * y*(q[2] - q[0] + y*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + y*(3.0*(q[1] - q[2]) + q[3] - q[0])));

			o_data[pos_idx] = value;
		}
	}



public:
	void bicubic_scalar(
			const PlaneData &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			PlaneData &o_data,					///< output values

			double i_shift_x = 0.0,				///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);


#if SWEET_USE_PLANE_SPECTRAL_SPACE
		o_data.physical_space_data_valid = true;
		o_data.spectral_space_data_valid = false;
#endif
	}


public:
	void bicubic_scalar(
			const PlaneData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,			///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,			///< y positions of interpolation points

			ScalarDataArray &o_data,				///< output values

			double i_shift_x = 0.0,            ///< shift in x for staggered grids
			double i_shift_y = 0.0				///< shift in y for staggered grids
	)
	{
		assert(res[0] > 0);
		assert(cached_scale_factor[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data, i_shift_x, i_shift_y);
	}


public:
	void bilinear_scalar(
			const PlaneData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			double *o_data,							///< output values
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

		std::size_t size = i_pos_x.number_of_elements;

		assert(size != 0);

		// iterate over all positions
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t pos_idx = 0; pos_idx < size; pos_idx++)
		{
			// load position to interpolate
			double pos_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*cached_scale_factor[0] + i_shift_x, (double)res[0]);
			double pos_y = wrapPeriodic(i_pos_y.scalar_data[pos_idx]*cached_scale_factor[1] + i_shift_y, (double)res[1]);

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

// TODO: Each 2nd wrapPeriodic is obsolete!!!
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
				p[0] = i_data.p_physical_get(idx_j, idx_i[0]);
				p[1] = i_data.p_physical_get(idx_j, idx_i[1]);

				q[kj] = p[0] + x*(p[1]-p[0]);

				idx_j = wrapPeriodic(idx_j+1, res[1]);
				//std::cout<< p[0] << p[1] << x << std::endl;
			}

			double value = q[0] + y*(q[1]-q[0]);

			o_data[pos_idx] = value;
		}
	}


public:
	void bilinear_scalar(
			const PlaneData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			ScalarDataArray &o_data,				///< output values
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data, i_shift_x, i_shift_y);
	}


public:
	void bilinear_scalar(
			const PlaneData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

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
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.planeDataConfig->physical_array_data_number_of_elements);

		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_shift_x, i_shift_y);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		o_data.physical_space_data_valid = true;
		o_data.spectral_space_data_valid = false;
#endif
	}


public:
	const ScalarDataArray bilinear_scalar(
			const PlaneData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		ScalarDataArray out(i_data.planeDataConfig->physical_array_data_number_of_elements);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}

public:
	const PlaneData bicubic_scalar(
			PlaneData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points
			//PlaneData* i_pos[2],	///< sampling position
			double i_shift_x = 0.0,
			double i_shift_y = 0.0
	)
	{
		PlaneData out(planeDataConfig);
		bicubic_scalar(i_data, i_pos_x, i_pos_y, out, i_shift_x, i_shift_y);
		return out;
	}
};




#endif /* SRC_INCLUDE_SWEET_PLANEDATASAMPLER_HPP_ */
