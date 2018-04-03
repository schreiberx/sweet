/*
 * sPHERESampler2D.hpp
 *
 *  Created on: 29 Mar 2018
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_SPHEREDATASAMPLER_HPP_
#define SRC_INCLUDE_SWEET_SPHEREDATASAMPLER_HPP_

#include <sweet/ScalarDataArray.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>


/**
 * this is a sampler class which provides method to provide
 * interpolated sampled values on 2D physical sphere data which
 * is provided by SphereData
 */
class SphereDataSampler
{
public:
	int res[2];						/// resolution of domain
	const SphereDataConfig *sphereDataConfig;

	std::vector<double> sampling_data;

private:
//	double cached_scale_factor[2];			/// cached parameters for sampling

	// number of extended latitudinal points (num_lat + 4)
	int ext_lat_M;

	// lookup table with latitudinal angles extended by 2 at front and back
	std::vector<double> phi_lookup;

	// distance between phi angles
	std::vector<double> phi_dist;


public:
	SphereDataSampler(
		SphereDataConfig *i_sphereDataConfig
	)
	{
		assert(i_sphereDataConfig != nullptr);

		sphereDataConfig = i_sphereDataConfig;
		setup(sphereDataConfig);
	}


	SphereDataSampler()
	{
		sphereDataConfig = nullptr;

		res[0] = -1;
		res[1] = -1;

//		cached_scale_factor[0] = -1;
//		cached_scale_factor[1] = -1;
	}



public:
	void setup(
		const SphereDataConfig *i_sphereDataConfig
	)
	{
		assert(i_sphereDataConfig != nullptr);
		sphereDataConfig = i_sphereDataConfig;

		res[0] = i_sphereDataConfig->physical_num_lon;
		res[1] = i_sphereDataConfig->physical_num_lat;

		ext_lat_M = sphereDataConfig->physical_num_lat+4;
		phi_lookup.resize(ext_lat_M);

		for (int i = 0; i < sphereDataConfig->physical_num_lat; i++)
			phi_lookup[i+2] = sphereDataConfig->lat[i];
		phi_lookup[1] = M_PI - phi_lookup[2];
		phi_lookup[0] = M_PI - phi_lookup[3];
		phi_lookup[ext_lat_M-2] = -M_PI - phi_lookup[ext_lat_M-3];
		phi_lookup[ext_lat_M-1] = -M_PI - phi_lookup[ext_lat_M-4];

	//	for (int i = 0; i < ext_lat_M; i++)
	//		std::cout << phi_lookup[i] << std::endl;

		phi_dist.resize(ext_lat_M-1);
		for (int i = 0; i < ext_lat_M-1; i++)
		{
			phi_dist[i] = phi_lookup[i] - phi_lookup[i+1];
			assert(phi_dist[i] > 0);
		}

	//	for (int i = 0; i < ext_lat_M-1; i++)
	//		std::cout << phi_dist[i] << std::endl;

	}



	void updateSamplingData(const SphereData &i_data)
	{
		i_data.request_data_physical();

		sampling_data.resize(sphereDataConfig->physical_num_lon*(sphereDataConfig->physical_num_lat+4));
//		for (int i = 0; i < sampling_data.size(); i++)
//			sampling_data[i] = -1;

		int num_lon = sphereDataConfig->physical_num_lon;
		int num_lat = sphereDataConfig->physical_num_lat;

		int num_lon_d2 = sphereDataConfig->physical_num_lon/2;

		assert((num_lon & 1) == 0);

#if 0
		((SphereData&)i_data).physical_update_lambda_array(
				[&](int x, int y, double &io_data)
				{
					io_data = x*10000+y;
				}
			);
		i_data.physical_print();
#endif

#if 1
		// first block
		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[i] = i_data.physical_space_data[num_lon + num_lon_d2 + i];
		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[num_lon_d2 + i] = i_data.physical_space_data[num_lon + i];

		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[num_lon + i] = i_data.physical_space_data[num_lon_d2 + i];
		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[num_lon + num_lon_d2 + i] = i_data.physical_space_data[i];
#endif

		for (int i = 0; i < num_lon*num_lat; i++)
			sampling_data[i + num_lon*2] = i_data.physical_space_data[i];


#if 1
		// last block
		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[num_lon*(num_lat+2) + i] = i_data.physical_space_data[num_lon*(num_lat-1) + num_lon_d2 + i];
		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[num_lon*(num_lat+2) + num_lon_d2 + i] = i_data.physical_space_data[num_lon*(num_lat-1) + i];

		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[num_lon*(num_lat+3) + i] = i_data.physical_space_data[num_lon*(num_lat-2) + num_lon_d2 + i];
		for (int i = 0; i < num_lon_d2; i++)
			sampling_data[num_lon*(num_lat+3) + num_lon_d2 + i] = i_data.physical_space_data[num_lon*(num_lat-2) + i];
#endif

#if 0
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (int j = 0; j < num_lat+4; j++)
		{
			for (int i = 0; i < num_lon; i++)
			{
				std::cout << sampling_data[i+j*num_lon] << "\t";
			}
			std::cout << std::endl;
		}
#endif
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
		/*
		 * TODO: replace this with efficient hardware operation (if available)
		 */
		while (i < 0)
			i += i_res;

		while (i >= i_res)
			i -= i_res;

		assert(i >= 0 && i < i_res);

		return i;
	}


	/*
	 * Compute determinant of 4x4 matrix
	 */
	static
	double get4x4Determinant(
			const double *i_m
	)
	{
		return	  i_m[0*4+3] * i_m[1*4+2] * i_m[2*4+1] * i_m[3*4+0] - i_m[0*4+2] * i_m[1*4+3] * i_m[2*4+1] * i_m[3*4+0] - i_m[0*4+3] * i_m[1*4+1] * i_m[2*4+2] * i_m[3*4+0] + i_m[0*4+1] * i_m[1*4+3] * i_m[2*4+2] * i_m[3*4+0]
				+ i_m[0*4+2] * i_m[1*4+1] * i_m[2*4+3] * i_m[3*4+0] - i_m[0*4+1] * i_m[1*4+2] * i_m[2*4+3] * i_m[3*4+0] - i_m[0*4+3] * i_m[1*4+2] * i_m[2*4+0] * i_m[3*4+1] + i_m[0*4+2] * i_m[1*4+3] * i_m[2*4+0] * i_m[3*4+1]
				+ i_m[0*4+3] * i_m[1*4+0] * i_m[2*4+2] * i_m[3*4+1] - i_m[0*4+0] * i_m[1*4+3] * i_m[2*4+2] * i_m[3*4+1] - i_m[0*4+2] * i_m[1*4+0] * i_m[2*4+3] * i_m[3*4+1] + i_m[0*4+0] * i_m[1*4+2] * i_m[2*4+3] * i_m[3*4+1]
				+ i_m[0*4+3] * i_m[1*4+1] * i_m[2*4+0] * i_m[3*4+2] - i_m[0*4+1] * i_m[1*4+3] * i_m[2*4+0] * i_m[3*4+2] - i_m[0*4+3] * i_m[1*4+0] * i_m[2*4+1] * i_m[3*4+2] + i_m[0*4+0] * i_m[1*4+3] * i_m[2*4+1] * i_m[3*4+2]
				+ i_m[0*4+1] * i_m[1*4+0] * i_m[2*4+3] * i_m[3*4+2] - i_m[0*4+0] * i_m[1*4+1] * i_m[2*4+3] * i_m[3*4+2] - i_m[0*4+2] * i_m[1*4+1] * i_m[2*4+0] * i_m[3*4+3] + i_m[0*4+1] * i_m[1*4+2] * i_m[2*4+0] * i_m[3*4+3]
				+ i_m[0*4+2] * i_m[1*4+0] * i_m[2*4+1] * i_m[3*4+3] - i_m[0*4+0] * i_m[1*4+2] * i_m[2*4+1] * i_m[3*4+3] - i_m[0*4+1] * i_m[1*4+0] * i_m[2*4+2] * i_m[3*4+3] + i_m[0*4+0] * i_m[1*4+1] * i_m[2*4+2] * i_m[3*4+3];
	}


	/*
	 * Solve 4x4 system of equation
	 */
	static
	void solve4x4SOE(
			const double *i_mat,	///< matrix
			const double *i_b,		///< rhs
			double *o_x				///< solution
	)
	{

		double minv[16] = {
			i_mat[1*4+2]*i_mat[2*4+3]*i_mat[3*4+1] - i_mat[1*4+3]*i_mat[2*4+2]*i_mat[3*4+1] + i_mat[1*4+3]*i_mat[2*4+1]*i_mat[3*4+2] - i_mat[1*4+1]*i_mat[2*4+3]*i_mat[3*4+2] - i_mat[1*4+2]*i_mat[2*4+1]*i_mat[3*4+3] + i_mat[1*4+1]*i_mat[2*4+2]*i_mat[3*4+3],
			i_mat[0*4+3]*i_mat[2*4+2]*i_mat[3*4+1] - i_mat[0*4+2]*i_mat[2*4+3]*i_mat[3*4+1] - i_mat[0*4+3]*i_mat[2*4+1]*i_mat[3*4+2] + i_mat[0*4+1]*i_mat[2*4+3]*i_mat[3*4+2] + i_mat[0*4+2]*i_mat[2*4+1]*i_mat[3*4+3] - i_mat[0*4+1]*i_mat[2*4+2]*i_mat[3*4+3],
			i_mat[0*4+2]*i_mat[1*4+3]*i_mat[3*4+1] - i_mat[0*4+3]*i_mat[1*4+2]*i_mat[3*4+1] + i_mat[0*4+3]*i_mat[1*4+1]*i_mat[3*4+2] - i_mat[0*4+1]*i_mat[1*4+3]*i_mat[3*4+2] - i_mat[0*4+2]*i_mat[1*4+1]*i_mat[3*4+3] + i_mat[0*4+1]*i_mat[1*4+2]*i_mat[3*4+3],
			i_mat[0*4+3]*i_mat[1*4+2]*i_mat[2*4+1] - i_mat[0*4+2]*i_mat[1*4+3]*i_mat[2*4+1] - i_mat[0*4+3]*i_mat[1*4+1]*i_mat[2*4+2] + i_mat[0*4+1]*i_mat[1*4+3]*i_mat[2*4+2] + i_mat[0*4+2]*i_mat[1*4+1]*i_mat[2*4+3] - i_mat[0*4+1]*i_mat[1*4+2]*i_mat[2*4+3],

			i_mat[1*4+3]*i_mat[2*4+2]*i_mat[3*4+0] - i_mat[1*4+2]*i_mat[2*4+3]*i_mat[3*4+0] - i_mat[1*4+3]*i_mat[2*4+0]*i_mat[3*4+2] + i_mat[1*4+0]*i_mat[2*4+3]*i_mat[3*4+2] + i_mat[1*4+2]*i_mat[2*4+0]*i_mat[3*4+3] - i_mat[1*4+0]*i_mat[2*4+2]*i_mat[3*4+3],
			i_mat[0*4+2]*i_mat[2*4+3]*i_mat[3*4+0] - i_mat[0*4+3]*i_mat[2*4+2]*i_mat[3*4+0] + i_mat[0*4+3]*i_mat[2*4+0]*i_mat[3*4+2] - i_mat[0*4+0]*i_mat[2*4+3]*i_mat[3*4+2] - i_mat[0*4+2]*i_mat[2*4+0]*i_mat[3*4+3] + i_mat[0*4+0]*i_mat[2*4+2]*i_mat[3*4+3],
			i_mat[0*4+3]*i_mat[1*4+2]*i_mat[3*4+0] - i_mat[0*4+2]*i_mat[1*4+3]*i_mat[3*4+0] - i_mat[0*4+3]*i_mat[1*4+0]*i_mat[3*4+2] + i_mat[0*4+0]*i_mat[1*4+3]*i_mat[3*4+2] + i_mat[0*4+2]*i_mat[1*4+0]*i_mat[3*4+3] - i_mat[0*4+0]*i_mat[1*4+2]*i_mat[3*4+3],
			i_mat[0*4+2]*i_mat[1*4+3]*i_mat[2*4+0] - i_mat[0*4+3]*i_mat[1*4+2]*i_mat[2*4+0] + i_mat[0*4+3]*i_mat[1*4+0]*i_mat[2*4+2] - i_mat[0*4+0]*i_mat[1*4+3]*i_mat[2*4+2] - i_mat[0*4+2]*i_mat[1*4+0]*i_mat[2*4+3] + i_mat[0*4+0]*i_mat[1*4+2]*i_mat[2*4+3],

			i_mat[1*4+1]*i_mat[2*4+3]*i_mat[3*4+0] - i_mat[1*4+3]*i_mat[2*4+1]*i_mat[3*4+0] + i_mat[1*4+3]*i_mat[2*4+0]*i_mat[3*4+1] - i_mat[1*4+0]*i_mat[2*4+3]*i_mat[3*4+1] - i_mat[1*4+1]*i_mat[2*4+0]*i_mat[3*4+3] + i_mat[1*4+0]*i_mat[2*4+1]*i_mat[3*4+3],
			i_mat[0*4+3]*i_mat[2*4+1]*i_mat[3*4+0] - i_mat[0*4+1]*i_mat[2*4+3]*i_mat[3*4+0] - i_mat[0*4+3]*i_mat[2*4+0]*i_mat[3*4+1] + i_mat[0*4+0]*i_mat[2*4+3]*i_mat[3*4+1] + i_mat[0*4+1]*i_mat[2*4+0]*i_mat[3*4+3] - i_mat[0*4+0]*i_mat[2*4+1]*i_mat[3*4+3],
			i_mat[0*4+1]*i_mat[1*4+3]*i_mat[3*4+0] - i_mat[0*4+3]*i_mat[1*4+1]*i_mat[3*4+0] + i_mat[0*4+3]*i_mat[1*4+0]*i_mat[3*4+1] - i_mat[0*4+0]*i_mat[1*4+3]*i_mat[3*4+1] - i_mat[0*4+1]*i_mat[1*4+0]*i_mat[3*4+3] + i_mat[0*4+0]*i_mat[1*4+1]*i_mat[3*4+3],
			i_mat[0*4+3]*i_mat[1*4+1]*i_mat[2*4+0] - i_mat[0*4+1]*i_mat[1*4+3]*i_mat[2*4+0] - i_mat[0*4+3]*i_mat[1*4+0]*i_mat[2*4+1] + i_mat[0*4+0]*i_mat[1*4+3]*i_mat[2*4+1] + i_mat[0*4+1]*i_mat[1*4+0]*i_mat[2*4+3] - i_mat[0*4+0]*i_mat[1*4+1]*i_mat[2*4+3],

			i_mat[1*4+2]*i_mat[2*4+1]*i_mat[3*4+0] - i_mat[1*4+1]*i_mat[2*4+2]*i_mat[3*4+0] - i_mat[1*4+2]*i_mat[2*4+0]*i_mat[3*4+1] + i_mat[1*4+0]*i_mat[2*4+2]*i_mat[3*4+1] + i_mat[1*4+1]*i_mat[2*4+0]*i_mat[3*4+2] - i_mat[1*4+0]*i_mat[2*4+1]*i_mat[3*4+2],
			i_mat[0*4+1]*i_mat[2*4+2]*i_mat[3*4+0] - i_mat[0*4+2]*i_mat[2*4+1]*i_mat[3*4+0] + i_mat[0*4+2]*i_mat[2*4+0]*i_mat[3*4+1] - i_mat[0*4+0]*i_mat[2*4+2]*i_mat[3*4+1] - i_mat[0*4+1]*i_mat[2*4+0]*i_mat[3*4+2] + i_mat[0*4+0]*i_mat[2*4+1]*i_mat[3*4+2],
			i_mat[0*4+2]*i_mat[1*4+1]*i_mat[3*4+0] - i_mat[0*4+1]*i_mat[1*4+2]*i_mat[3*4+0] - i_mat[0*4+2]*i_mat[1*4+0]*i_mat[3*4+1] + i_mat[0*4+0]*i_mat[1*4+2]*i_mat[3*4+1] + i_mat[0*4+1]*i_mat[1*4+0]*i_mat[3*4+2] - i_mat[0*4+0]*i_mat[1*4+1]*i_mat[3*4+2],
			i_mat[0*4+1]*i_mat[1*4+2]*i_mat[2*4+0] - i_mat[0*4+2]*i_mat[1*4+1]*i_mat[2*4+0] + i_mat[0*4+2]*i_mat[1*4+0]*i_mat[2*4+1] - i_mat[0*4+0]*i_mat[1*4+2]*i_mat[2*4+1] - i_mat[0*4+1]*i_mat[1*4+0]*i_mat[2*4+2] + i_mat[0*4+0]*i_mat[1*4+1]*i_mat[2*4+2]
		};

		for (int j = 0; j < 4; j++)
		{
			o_x[j] = 0;
			for (int i = 0; i < 4; i++)
				o_x[j] += minv[j*4+i]*i_b[i];
		}


		double inv_det = 1.0/get4x4Determinant(i_mat);
		for (int i = 0; i < 4; i++)
			o_x[i] *= inv_det;
	}


public:
	void bicubic_scalar(
			const SphereData &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data						///< output values
	)
	{
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		updateSamplingData(i_data);

		int num_lon = sphereDataConfig->physical_num_lon;
		int num_lat = sphereDataConfig->physical_num_lat;

		double s_lon = (double)i_data.sphereDataConfig->physical_num_lon / (2.0*M_PI);

		double L = -(-M_PI*0.5 - M_PI/ext_lat_M*1.5);
		// total size of lat field (M_PI + extension)
		// divided by number of cells
		//double s = (M_PI+M_PI/ext_lat_M*3)/(double)(ext_lat_M-1);
		double inv_s = (double)(ext_lat_M-1)/(M_PI+M_PI/ext_lat_M*3);

		// iterate over all positions in parallel
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t pos_idx = 0; pos_idx < i_pos_x.number_of_elements; pos_idx++)
		{
			// compute array position
			double array_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*s_lon, (double)res[0]);
			// compute position relative in cell \in [0;1]
			double cell_x = array_x - std::floor(array_x);
			assert(cell_x >= 0);
			assert(cell_x <= 1);

			// compute array index
			int array_idx_x = std::floor(array_x);
			assert(array_idx_x >= 0);
			assert(array_idx_x < sphereDataConfig->physical_num_lon);

			// estimate array index for latitude
			double phi = i_pos_y.scalar_data[pos_idx];
			int est_lat_idx = (L - phi)*inv_s;

			assert(est_lat_idx >= 0);
			assert(est_lat_idx < ext_lat_M-1);

			if (phi_lookup[est_lat_idx] < phi)
				est_lat_idx--;
			else if (phi_lookup[est_lat_idx+1] > phi)
				est_lat_idx++;

			int array_idx_y = est_lat_idx;
			assert(array_idx_y >= 0);
			assert(array_idx_y < ext_lat_M);

			double cell_y = (phi - phi_lookup[array_idx_y+1]) / phi_dist[array_idx_y];
//			std::cout << cell_y << std::endl;
			// flip since the coordinate system is also flipped!
			cell_y = 1.0-cell_y;
//			std::cout << cell_y << std::endl;
			assert(cell_y >= 0);
			assert(cell_y <= 1);

			assert(phi_lookup[array_idx_y] >= phi);
			assert(phi_lookup[array_idx_y+1] <= phi);


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

				q[kj] = p[1] + 0.5 * cell_x*(p[2] - p[0] + cell_x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + cell_x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

				idx_j++;
			}

			//phi_dist[array_idx_y]
			double mat[16];
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					//double x = j;//phi_dist[array_idx_y-1+j];
					double x = phi_lookup[array_idx_y-1+j];
					mat[j*4+i] = std::pow(x, i);
				}
			}

			double x[4];
			solve4x4SOE(mat, q, x);

#if 1
			//double y = cell_y+1.0;
			double y = phi;
			//y = phi_dist[array_idx_y-1]+y*;
			double value = x[0] + y*(x[1] + y*(x[2] + y*x[3]));

#else
			double value = q[1] + 0.5 * cell_y*(q[2] - q[0] + cell_y*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + cell_y*(3.0*(q[1] - q[2]) + q[3] - q[0])));
#endif

			o_data[pos_idx] = value;
		}
	}


public:
	void bicubic_scalar(
			const SphereData &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			SphereData &o_data					///< output values
	)
	{
		assert(i_data.sphereDataConfig->physical_array_data_number_of_elements == o_data.sphereDataConfig->physical_array_data_number_of_elements);
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == (std::size_t)o_data.sphereDataConfig->physical_array_data_number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data);

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		o_data.physical_space_data_valid = true;
		o_data.spectral_space_data_valid = false;
#endif
	}


public:
	void bicubic_scalar(
			const SphereData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,			///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,			///< y positions of interpolation points

			ScalarDataArray &o_data				///< output values
	)
	{
//		assert ((std::size_t)i_data.sphereDataConfig->physical_array_data_number_of_elements == (std::size_t)o_data.number_of_elements);
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == o_data.number_of_elements);

		bicubic_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data);
	}



public:
	void bilinear_scalar(
			const SphereData &i_data,			///< sampling data

			const ScalarDataArray &i_pos_x,		///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,		///< y positions of interpolation points

			double *o_data						///< output values
	)
	{
		assert(res[0] > 0);
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);

		updateSamplingData(i_data);

		int num_lon = sphereDataConfig->physical_num_lon;
		int num_lat = sphereDataConfig->physical_num_lat;


		double s_lon = (double)i_data.sphereDataConfig->physical_num_lon / (2.0*M_PI);

		double L = -(-M_PI*0.5 - M_PI/ext_lat_M*1.5);
		// total size of lat field (M_PI + extension)
		// divided by number of cells
		//double s = (M_PI+M_PI/ext_lat_M*3)/(double)(ext_lat_M-1);
		double inv_s = (double)(ext_lat_M-1)/(M_PI+M_PI/ext_lat_M*3);

		// iterate over all positions in parallel
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t pos_idx = 0; pos_idx < i_pos_x.number_of_elements; pos_idx++)
		{
			// compute array position
			double array_x = wrapPeriodic(i_pos_x.scalar_data[pos_idx]*s_lon, (double)res[0]);
			// compute position relative in cell \in [0;1]
			double cell_x = array_x - std::floor(array_x);
			assert(cell_x >= 0);
			assert(cell_x <= 1);

			// compute array index
			int array_idx_x = std::floor(array_x);
			assert(array_idx_x >= 0);
			assert(array_idx_x < sphereDataConfig->physical_num_lon);

			// estimate array index for latitude
			double phi = i_pos_y.scalar_data[pos_idx];
			int est_lat_idx = (L - phi)*inv_s;

			assert(est_lat_idx >= 0);
			assert(est_lat_idx < ext_lat_M-1);

			if (phi_lookup[est_lat_idx] < phi)
				est_lat_idx--;
			else if (phi_lookup[est_lat_idx+1] > phi)
				est_lat_idx++;

			int array_idx_y = est_lat_idx;
			assert(array_idx_y >= 0);
			assert(array_idx_y < ext_lat_M);

			double cell_y = (phi - phi_lookup[array_idx_y+1]) / phi_dist[array_idx_y];
			// flip since the coordinate system is also flipped!
			cell_y = 1.0-cell_y;
			assert(cell_y >= 0);
			assert(cell_y <= 1);

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

				q[kj] = p[0] + cell_x*(p[1]-p[0]);

				idx_j++;
			}
			double value = q[0] + cell_y*(q[1]-q[0]);

			o_data[pos_idx] = value;
		}
	}

public:
	void bilinear_scalar(
			const SphereData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			ScalarDataArray &o_data				///< output values
	)
	{
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == (std::size_t)o_data.number_of_elements);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.scalar_data);
	}


public:
	void bilinear_scalar(
			const SphereData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y,				///< y positions of interpolation points

			SphereData &o_data				///< output values
	)
	{
		assert(i_pos_x.number_of_elements == i_pos_y.number_of_elements);
		assert(i_pos_x.number_of_elements == (std::size_t)o_data.sphereDataConfig->physical_array_data_number_of_elements);

		bilinear_scalar(i_data, i_pos_x, i_pos_y, o_data.physical_space_data);

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		o_data.physical_space_data_valid = true;
		o_data.spectral_space_data_valid = false;
#endif
	}


public:
	const ScalarDataArray bilinear_scalar(
			const SphereData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y				///< y positions of interpolation points
	)
	{
		ScalarDataArray out(i_data.sphereDataConfig->physical_array_data_number_of_elements);
		bilinear_scalar(i_data, i_pos_x, i_pos_y, out);
		return out;
	}

public:
	const SphereData bicubic_scalar(
			SphereData &i_data,				///< sampling data

			const ScalarDataArray &i_pos_x,				///< x positions of interpolation points
			const ScalarDataArray &i_pos_y				///< y positions of interpolation points
	)
	{
		SphereData out(sphereDataConfig);
		bicubic_scalar(i_data, i_pos_x, i_pos_y, out);
		return out;
	}
};




#endif /* SRC_INCLUDE_SWEET_SPHEREDATASAMPLER_HPP_ */
