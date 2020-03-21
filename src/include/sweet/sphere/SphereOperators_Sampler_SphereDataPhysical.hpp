/*
 * SphereDataSampler2D.hpp
 *
 *  Created on: 29 Mar 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHEREDATASAMPLER_HPP_
#define SRC_INCLUDE_SWEET_SPHEREDATASAMPLER_HPP_

#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/ScalarDataArray.hpp>
#include <libmath/interpolation.hpp>


/**
 * this is a sampler class which provides method to provide
 * interpolated sampled values on 2D physical sphere data which
 * is provided by SphereDataPhysical
 */
class SphereOperators_Sampler_SphereDataPhysical
{
public:
	int res[2];						/// resolution of domain
	const SphereData_Config *sphereDataConfig;

	std::vector<double> sampling_data;

private:
//	double cached_scale_factor[2];			/// cached parameters for sampling

	// number of extended latitudinal points (num_lat + 4)
	int ext_lat_M;

	// lookup table with latitudinal angles extended by 2 at front and back
	std::vector<double> phi_lookup;

#if 0
	// distance between phi angles
	std::vector<double> phi_dist;

	// storage for inverse matrices
	std::vector<double> inv_matrices;
#endif

public:
	SphereOperators_Sampler_SphereDataPhysical(
		SphereData_Config *i_sphereDataConfig
	)
	{
		assert(i_sphereDataConfig != nullptr);

		sphereDataConfig = i_sphereDataConfig;
		setup(sphereDataConfig);
	}


	SphereOperators_Sampler_SphereDataPhysical()
	{
		sphereDataConfig = nullptr;

		res[0] = -1;
		res[1] = -1;
	}



public:
	void setup(
		const SphereData_Config *i_sphereDataConfig
	)
	{
		assert(i_sphereDataConfig != nullptr);
		sphereDataConfig = i_sphereDataConfig;

		res[0] = i_sphereDataConfig->physical_num_lon;
		res[1] = i_sphereDataConfig->physical_num_lat;


		/*
		 * Use extended lat/lon lookup table to avoid if brances
		 */
		ext_lat_M = sphereDataConfig->physical_num_lat+4;
		phi_lookup.resize(ext_lat_M);

		for (int i = 0; i < sphereDataConfig->physical_num_lat; i++)
			phi_lookup[i+2] = sphereDataConfig->lat[i];

		phi_lookup[0] = M_PI - phi_lookup[3];
		phi_lookup[1] = M_PI - phi_lookup[2];
		phi_lookup[ext_lat_M-1] = -M_PI - phi_lookup[ext_lat_M-4];
		phi_lookup[ext_lat_M-2] = -M_PI - phi_lookup[ext_lat_M-3];

#if 0
		for (int i = 0; i < sphereDataConfig->physical_num_lat+4; i++)
			std::cout << phi_lookup[i] << std::endl;
		exit(-1);
#endif


		/*
		 * Test for proper cubic interpolation
		 */
		double a[4] = {1.0, 2.0, 3.0, 4.0};

		double x0 = 1.3;

		auto f = [&](double x) -> double
		{
			return a[0] + x*a[1] + x*x*a[2] + x*x*x*a[3];
		};

		double x[4] = {0.1, 0.2, 0.4, 0.8};
		double y[4] = {f(x[0]), f(x[1]), f(x[2]), f(x[3])};

		double retval = interpolation_hermite_nonequidistant<4>(x, y, x0);

		if (std::abs(retval - f(x0)) > 1e-10)
			FatalError("Cubic interpolation Buggy!!!");
	}



	void updateSamplingData(
			const SphereData_Physical &i_data,
			bool i_velocity_sampling = false
	)
	{
		sampling_data.resize(sphereDataConfig->physical_num_lon*(sphereDataConfig->physical_num_lat+4));

		int num_lon = sphereDataConfig->physical_num_lon;
		int num_lat = sphereDataConfig->physical_num_lat;

		int num_lon_d2 = sphereDataConfig->physical_num_lon/2;

		assert((num_lon & 1) == 0);

		if (!i_velocity_sampling)
		{
			// first block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[i] = i_data.physical_space_data[num_lon + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon_d2 + i] = i_data.physical_space_data[num_lon + i];

			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon + i] = i_data.physical_space_data[num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon + num_lon_d2 + i] = i_data.physical_space_data[i];
		}
		else
		{
			// first block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[i] = -i_data.physical_space_data[num_lon + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon_d2 + i] = -i_data.physical_space_data[num_lon + i];

			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon + i] = -i_data.physical_space_data[num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon + num_lon_d2 + i] = -i_data.physical_space_data[i];
		}

		for (int i = 0; i < num_lon*num_lat; i++)
			sampling_data[i + num_lon*2] = i_data.physical_space_data[i];


		if (!i_velocity_sampling)
		{
			// last block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+2) + i] = i_data.physical_space_data[num_lon*(num_lat-1) + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+2) + num_lon_d2 + i] = i_data.physical_space_data[num_lon*(num_lat-1) + i];

			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+3) + i] = i_data.physical_space_data[num_lon*(num_lat-2) + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+3) + num_lon_d2 + i] = i_data.physical_space_data[num_lon*(num_lat-2) + i];
		}
		else
		{
			// last block
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+2) + i] = -i_data.physical_space_data[num_lon*(num_lat-1) + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+2) + num_lon_d2 + i] = -i_data.physical_space_data[num_lon*(num_lat-1) + i];

			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+3) + i] = -i_data.physical_space_data[num_lon*(num_lat-2) + num_lon_d2 + i];
			for (int i = 0; i < num_lon_d2; i++)
				sampling_data[num_lon*(num_lat+3) + num_lon_d2 + i] = -i_data.physical_space_data[num_lon*(num_lat-2) + i];
		}
	}



public:
	/**
	 * wrap the position i in a periodic domain of size i_res
	 */
	template <typename T>
	static
	inline
	T wrapPeriodic(T i, T i_res)
	{
#if 1

		i = fmod(i, i_res);
		if (i < 0)
			i += i_res;

#else

		while (i < 0)
			i += i_res;

		while (i >= i_res)
			i -= i_res;
#endif

		assert(i >= 0 && i < i_res);

		return i;
	}


public:
	void bicubic_scalar(
			const SphereData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values
			bool i_velocity_sampling,
			bool i_limiter			///< Use limiter for interpolation to avoid unphysical local extrema
	)
	{
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		updateSamplingData(i_data, i_velocity_sampling);

		int num_lon = sphereDataConfig->physical_num_lon;
		int num_lat = sphereDataConfig->physical_num_lat;

		double s_lon = (double)i_data.sphereDataConfig->physical_num_lon / (2.0*M_PI);

		double L = -(-M_PI*0.5 - M_PI/ext_lat_M*1.5);
		// total size of lat field (M_PI + extension)
		// divided by number of cells
		//double s = (M_PI+M_PI/ext_lat_M*3)/(double)(ext_lat_M-1);
		double inv_s = (double)(ext_lat_M-1)/(M_PI+M_PI/ext_lat_M*3);

		// iterate over all positions in parallel
#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (std::size_t pos_idx = 0; pos_idx < i_pos_x.number_of_elements; pos_idx++)
		{
			/*
			 * Compute X array position
			 */

			double array_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*s_lon, (double)res[0]);
			// compute position relative in cell \in [0;1]
			double cell_rel_x = array_x - std::floor(array_x);
			assert(cell_rel_x >= 0);
			assert(cell_rel_x <= 1);

			// compute array index
			int array_idx_x = std::floor(array_x);
			assert(array_idx_x >= 0);
			assert(array_idx_x < sphereDataConfig->physical_num_lon);

			/*
			 * Compute Y array position
			 */
			// estimate array index for latitude
			double phi = i_pos_y.scalar_data[pos_idx];
			int est_lat_idx = (L - phi)*inv_s;

#if SWEET_DEBUG
			if (!(est_lat_idx >= 0))
			{
				std::cout << "est_lat_idx: " << est_lat_idx << std::endl;
				std::cout << "L: " << L << std::endl;
				std::cout << "phi: " << phi << std::endl;
				std::cout << "inv_s: " << inv_s << std::endl;
				FatalError("est_lat_idx");
			}
#endif
			assert(est_lat_idx >= 1);
			assert(est_lat_idx < ext_lat_M-1);

			if (phi_lookup[est_lat_idx] < phi)
				est_lat_idx--;
			else if (phi_lookup[est_lat_idx+1] > phi)
				est_lat_idx++;

			int array_idx_y = est_lat_idx;
			assert(array_idx_y >= 0);
			assert(array_idx_y < ext_lat_M);


			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */

			// precompute x-position indices since they are reused 4 times
			int idx_i[4];
			{
				idx_i[0] = wrapPeriodic(array_idx_x-1, res[0]);
				idx_i[1] = wrapPeriodic(array_idx_x+0, res[0]);
				idx_i[2] = wrapPeriodic(array_idx_x+1, res[0]);
				idx_i[3] = wrapPeriodic(array_idx_x+2, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = array_idx_y-1;

			double q[4];
			for (int kj = 0; kj < 4; kj++)
			{
				assert(idx_j >= 0);
				assert(idx_j < num_lat+4);
				double p[4];

				p[0] = sampling_data[idx_j*num_lon + idx_i[0]];
				p[1] = sampling_data[idx_j*num_lon + idx_i[1]];
				p[2] = sampling_data[idx_j*num_lon + idx_i[2]];
				p[3] = sampling_data[idx_j*num_lon + idx_i[3]];

				q[kj] = interpolation_hermite_equidistant<4>(p, cell_rel_x+1.0);

				if (i_limiter)
				{
					double max = std::max(p[1], p[2]);
					double min = std::min(p[1], p[2]);

					q[kj] = std::min(q[kj], max);
					q[kj] = std::max(q[kj], min);
				}

				idx_j++;
			}

			double value = interpolation_hermite_nonequidistant<4>(&phi_lookup[array_idx_y-1], q, phi);

			if (i_limiter)
			{
				double max = std::max(q[1], q[2]);
				double min = std::min(q[1], q[2]);

				value = std::min(value, max);
				value = std::max(value, min);
			}

			o_data[pos_idx] = value;
		}
	}


public:
	void bicubic_scalar(
			const SphereData_Physical &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			SphereData_Physical &o_data,					///< output values
			bool i_velocity_sampling,
			bool i_limiter			///< Use limiter for interpolation to avoid unphysical local extrema
	)
	{
		assert(i_data.sphereDataConfig->physical_array_data_number_of_elements == o_data.sphereDataConfig->physical_array_data_number_of_elements);
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == (std::size_t)o_data.sphereDataConfig->physical_array_data_number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data,  i_velocity_sampling, i_limiter);
	}


public:
	void bicubic_scalar(
			const SphereData_Physical &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,			///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,			///< y positions of interpolation points

			ScalarDataArray &o_data,				///< output values
			bool i_velocity_sampling,				///< swap sign for velocities,
			bool i_limiter			///< Use limiter for interpolation to avoid unphysical local extrema
	)
	{
//		assert ((std::size_t)i_data.sphereDataConfig->physical_array_data_number_of_elements == (std::size_t)o_data.number_of_elements);
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data, i_velocity_sampling, i_limiter);
	}



public:
	const ScalarDataArray bicubic_scalar(
			const SphereData_Physical &i_data,		///< sampling data

			const ScalarDataArray &i_pos_x,			///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,			///< y positions of interpolation points

			bool i_velocity_sampling,				///< swap sign for velocities
			bool i_limiter			///< Use limiter for interpolation to avoid unphysical local extrema
	)
	{
		ScalarDataArray out(i_data.sphereDataConfig->physical_array_data_number_of_elements);
		bicubic_scalar(i_data, i_pos_x, i_pos_y, out, i_velocity_sampling, i_limiter);
		return out;
	}

public:
	void bilinear_scalar(
			const SphereData_Physical &i_data,	///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data,						///< output values
			bool i_velocity_sampling			///< swap sign for velocities,
	)
	{
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		// copy the data to an internal buffer including halo layers
		updateSamplingData(i_data, i_velocity_sampling);

		int num_lon = sphereDataConfig->physical_num_lon;
		int num_lat = sphereDataConfig->physical_num_lat;

		// longitude spacing
		double s_lon = (double)i_data.sphereDataConfig->physical_num_lon / (2.0*M_PI);

		double L = -(-M_PI*0.5 - M_PI/ext_lat_M*1.5);
		// total size of lat field (M_PI + extension)
		// divided by number of cells
		//double s = (M_PI+M_PI/ext_lat_M*3)/(double)(ext_lat_M-1);
		double inv_s = (double)(ext_lat_M-1)/(M_PI+M_PI/ext_lat_M*3);

		// iterate over all positions in parallel
		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t pos_idx = 0; pos_idx < i_pos_x.number_of_elements; pos_idx++)
		{
			/*
			 * Compute X information
			 */
			double array_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*s_lon, (double)res[0]);

			// compute position relative in cell \in [0;1]
			double cell_rel_x = array_x - std::floor(array_x);
			assert(cell_rel_x >= 0);
			assert(cell_rel_x <= 1);

			// compute array index
			int array_idx_x = std::floor(array_x);
			assert(array_idx_x >= 0);
			assert(array_idx_x < sphereDataConfig->physical_num_lon);

			/*
			 * Compute Y information
			 *
			 * This is done via a lookup into phi_lookup since these
			 * coordinates are not equidistantly spaced, but close to it
			 */
			// estimate array index for latitude
			double phi = i_pos_y.scalar_data[pos_idx];
			int est_lat_idx = (L - phi)*inv_s;
#if SWEET_DEBUG
			if (!(est_lat_idx >= 0))
			{
				std::cout << "pos_idx: " << pos_idx << std::endl;
				std::cout << "est_lat_idx: " << est_lat_idx << std::endl;
				std::cout << "L: " << L << std::endl;
				std::cout << "phi: " << phi << std::endl;
				std::cout << "inv_s: " << inv_s << std::endl;
				FatalError("est_lat_idx");
			}
#endif
			assert(est_lat_idx >= 0);
			assert(est_lat_idx < ext_lat_M-1);

			if (phi_lookup[est_lat_idx] < phi)
				est_lat_idx--;
			else if (phi_lookup[est_lat_idx+1] > phi)
				est_lat_idx++;

			int array_idx_y = est_lat_idx;
			assert(array_idx_y >= 0);
			assert(array_idx_y < ext_lat_M);
/*
			// compute relative position in cell
			double cell_rel_y = (phi - phi_lookup[array_idx_y+1]) / phi_dist[array_idx_y];
			// flip since the coordinate system is also flipped!
			cell_rel_y = 1.0-cell_rel_y;
			assert(cell_rel_y >= 0);
			assert(cell_rel_y <= 1);
*/
			assert(phi_lookup[array_idx_y] >= phi);
			assert(phi_lookup[array_idx_y+1] <= phi);


			/**
			 * See http://www.paulinternet.nl/?page=bicubic
			 */

			// precompute x-position indices since they are reused 4 times
			int idx_i[2];
			{
				idx_i[0] = wrapPeriodic(array_idx_x, res[0]);
				idx_i[1] = wrapPeriodic(array_idx_x+1, res[0]);
			}

			/**
			 * iterate over rows and interpolate over the columns in the x direction
			 */
			// start at this row
			int idx_j = array_idx_y;

			double q[2];
			for (int kj = 0; kj < 2; kj++)
			{
				assert(idx_j >= 0);
				assert(idx_j < num_lat+4);
				double p[2];

				p[0] = sampling_data[idx_j*num_lon + idx_i[0]];
				p[1] = sampling_data[idx_j*num_lon + idx_i[1]];

				q[kj] = interpolation_hermite_equidistant<2>(p, cell_rel_x);

				idx_j++;
			}

			// interpolation in y direction
			double value = interpolation_hermite_nonequidistant<2>(&phi_lookup[array_idx_y], q, phi);

			o_data[pos_idx] = value;
		}
	}

public:
	void bilinear_scalar(
			const SphereData_Physical &i_data,	///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			ScalarDataArray &o_data,			///< output values
			bool i_velocity_sampling			///< swap sign for velocities
	)
	{
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == (std::size_t)o_data.number_of_elements);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data, i_velocity_sampling);
	}


public:
	void bilinear_scalar(
			const SphereData_Physical &i_data,	///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			SphereData_Physical &o_data,		///< output values
			bool i_velocity_sampling			///< swap sign for velocities
	)
	{
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == (std::size_t)o_data.sphereDataConfig->physical_array_data_number_of_elements);

		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data, i_velocity_sampling);
	}


public:
	const ScalarDataArray bilinear_scalar(
			const SphereData_Physical &i_data,	///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			bool i_velocity_sampling			///< swap sign for velocities
	)
	{
		ScalarDataArray out(i_data.sphereDataConfig->physical_array_data_number_of_elements);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, out, i_velocity_sampling);
		return out;
	}


};




#endif /* SRC_INCLUDE_SWEET_SPHEREDATASAMPLER_HPP_ */
