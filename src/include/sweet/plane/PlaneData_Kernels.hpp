/*
 * PlaneData_Kernels.hpp
 *
 *  Created on: 19 Oct 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_PLANEDATA_KERNELS_HPP_
#define SRC_INCLUDE_SWEET_PLANE_PLANEDATA_KERNELS_HPP_

namespace sweet
{

class PlaneData_Kernels
{

//#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	int kernel_size = -1;
	double *kernel_data = nullptr;
	int kernel_id = -1;
//#endif




public:
	template <int S>
	void kernel_stencil_setup(
			const double i_kernel_array[S][S],
			double i_scale,

			const PlaneDataConfig *planeDataConfig,
			double *o_physical_data
	)
	{
//#if !SWEET_USE_PLANE_SPECTRAL_SPACE

		kernel_size = S;
		kernel_data = MemBlockAlloc::alloc<double>(sizeof(double)*S*S);
		for (int y = 0; y < S; y++)
			for (int x = 0; x < S; x++)
				kernel_data[y*S+x] = i_kernel_array[S-1-y][x];

		for (int i = 0; i < S*S; i++)
			kernel_data[i] *= i_scale;


		if (S == 3)
		{
			kernel_id = PlaneData_Kernels::get_kernel_mask3x3(
					kernel_data[0] != 0,
					kernel_data[1] != 0,
					kernel_data[2] != 0,
					kernel_data[3] != 0,
					kernel_data[4] != 0,
					kernel_data[5] != 0,
					kernel_data[6] != 0,
					kernel_data[7] != 0,
					kernel_data[8] != 0
				);
		}
		else
		{
			kernel_id = -1;
		}

//#else

#define KERNEL_physical_set(j, i, value)		\
	o_physical_data[(j)*planeDataConfig->physical_data_size[0]+(i)] = (value);


		double inv_kernel_array[S][S];

		for (int j = 0; j < S; j++)
			for (int i = 0; i < S; i++)
				inv_kernel_array[j][i] = i_kernel_array[j][S-i-1]*i_scale;

		// assure symmetric kernel
		assert((S & 1) == 1);

		// radius of kernel (half size)
		std::size_t R = S>>1;

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < planeDataConfig->physical_array_data_number_of_elements; i++)
			o_physical_data[i] = 0;

		// left lower corner
		//kernel_cart[    0:R[0]+1,       0:R[0]+1    ] = conv_kernel[    R[0]:,  R[0]:   ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				int cy = ky-R;
				int cx = kx-R;

				if (0 > cx || (int)planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || (int)planeDataConfig->physical_data_size[1] <= cy)
					continue;

				KERNEL_physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// right bottom corner
		//kernel_cart[    0:R[0]+1,       res-R[0]:   ] = conv_kernel[    R[0]:,  0:R[0]  ]
		for (std::size_t ky = R; ky < S; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				int cy = ky-R;
				int cx = planeDataConfig->physical_data_size[0] - R + kx;

				if (0 > cx || (int)planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || (int)planeDataConfig->physical_data_size[1] <= cy)
					continue;

				KERNEL_physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// left top corner
		//kernel_cart[    res-R[0]:,      0:R[0]+1   ] = conv_kernel[    0:R[0],  R[0]:  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = R; kx < S; kx++)
			{
				// coordinates in Cartesian space
				int cy = planeDataConfig->physical_data_size[1] - R + ky;
				int cx = kx-R;

				if (0 > cx || (int)planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || (int)planeDataConfig->physical_data_size[1] <= cy)
					continue;

				KERNEL_physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}


		// right top corner
		//kernel_cart[    res-R[0]:,      res-R[0]:   ] = conv_kernel[    0:R[0], 0:R[0]  ]
		for (std::size_t ky = 0; ky < R; ky++)
			for (std::size_t kx = 0; kx < R; kx++)
			{
				// coordinates in Cartesian space
				int cy = planeDataConfig->physical_data_size[1] - R + ky;
				int cx = planeDataConfig->physical_data_size[0] - R + kx;

				if (0 > cx || (int)planeDataConfig->physical_data_size[0] <= cx)
					continue;
				if (0 > cy || (int)planeDataConfig->physical_data_size[1] <= cy)
					continue;

				KERNEL_physical_set(cy, cx, inv_kernel_array[ky][kx]);
			}


#undef KERNEL_physical_set

//#endif
	}

public:
	constexpr
	static
	int get_kernel_mask3x3(
			int i_0,
			int i_1,
			int i_2,
			int i_3,
			int i_4,
			int i_5,
			int i_6,
			int i_7,
			int i_8
	)
	{
		return
				(i_0 << 0) |
				(i_1 << 1) |
				(i_2 << 2) |
				(i_3 << 3) |
				(i_4 << 4) |
				(i_5 << 5) |
				(i_6 << 6) |
				(i_7 << 7) |
				(i_8 << 8);
	}


///#if !SWEET_USE_PLANE_SPECTRAL_SPACE

public:
	void kernel_apply(
			int res_x,
			int res_y,
			double *i_data,

			double *o_data
	)	const
	{

		/**
		 * TODO: optimize this!!!
		 *
		 *  - cache blocking
		 *  - if branching elimination
		 *  - etc.....
		 */

		if (kernel_size == 3)
		{
			switch (kernel_id)
			{
			case get_kernel_mask3x3(0, 0, 0, 1, 0, 1, 0, 0, 0):	// (X, 0, X)

				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (int y = 0; y < res_y; y++)
				{
					for (int x = 0; x < res_x; x++)
					{
						double &data_out = o_data[y*res_x+x];
						data_out = 0;

						int pos_y = y;
						assert(pos_y >= 0 && pos_y < res_y);

						if (x > 0 && x < res_x-1)
						{
							double *kernel_scalar_ptr = &kernel_data[3];
							double *data_scalar_ptr = &i_data[pos_y*res_x+x-1];

							data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
							data_out += kernel_scalar_ptr[2]*data_scalar_ptr[2];
						}
						else
						{
							for (int i = -1; i <= 1; i+=2)
							{
								int pos_x = x+i;
								pos_x -= (pos_x >= res_x ? res_x : 0);
								pos_x += (pos_x < 0 ? res_x : 0);
								int idx = i+4;
								double kernel_scalar = kernel_data[idx];
								double data_scalar = i_data[pos_y*res_x+pos_x];

								data_out += kernel_scalar*data_scalar;
							}
						}
					}
				}
				break;



			case get_kernel_mask3x3(0, 0, 0, 1, 1, 1, 0, 0, 0):	// (X, X, X)
					SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
					for (int y = 0; y < res_y; y++)
					{
						for (int x = 0; x < res_x; x++)
						{
							double &data_out = o_data[y*res_x+x];
							data_out = 0;

							int pos_y = y+res_y;
							pos_y -= (pos_y >= res_y ? res_y : 0);
							pos_y += (pos_y < 0 ? res_y : 0);

							assert(pos_y >= 0 && pos_y < res_y);

							if (x > 0 && x < res_x-1)
							{
								double *kernel_scalar_ptr = &kernel_data[3];
								double *data_scalar_ptr = &i_data[pos_y*res_x+x-1];

								data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
								data_out += kernel_scalar_ptr[1]*data_scalar_ptr[1];
								data_out += kernel_scalar_ptr[2]*data_scalar_ptr[2];
							}
							else
							{
								for (int i = -1; i <= 1; i++)
								{
									int pos_x = (x+i);
									pos_x -= (pos_x >= res_x ? res_x : 0);
									pos_x += (pos_x < 0 ? res_x : 0);
									int idx = i+4;
									double kernel_scalar = kernel_data[idx];
									double data_scalar = i_data[pos_y*res_x+pos_x];

									data_out += kernel_scalar*data_scalar;
								}
							}
						}
					}
					break;

			case get_kernel_mask3x3(0, 1, 0, 0, 0, 0, 0, 1, 0):	// (X, 0, X)^T

					SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
					for (int y = 0; y < res_y; y++)
					{
						for (int x = 0; x < res_x; x++)
						{
							double &data_out = o_data[y*res_x+x];
							data_out = 0;

							if (y > 0 && y < res_y-1)
							{
								double *kernel_scalar_ptr = &kernel_data[1];
								double *data_scalar_ptr = &i_data[(y-1)*res_x+x];

								data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
								data_out += kernel_scalar_ptr[6]*data_scalar_ptr[2*res_x];
							}
							else
							{
								for (int j = -1; j <= 1; j+=2)
								{
									int pos_y = y+j;

									pos_y -= (pos_y >= res_y ? res_y : 0);
									pos_y += (pos_y < 0 ? res_y : 0);

									int idx = (j+1)*3+1;

									double kernel_scalar = kernel_data[idx];
									double data_scalar = i_data[pos_y*res_x+x];

									data_out += kernel_scalar*data_scalar;
								}
							}
						}
					}
					break;

			case get_kernel_mask3x3(0, 1, 0, 0, 1, 0, 0, 1, 0):	// (X, 0, X)^T

					SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
					for (int y = 0; y < res_y; y++)
					{
						for (int x = 0; x < res_x; x++)
						{
							double &data_out = o_data[y*res_x+x];
							data_out = 0;

							if (y > 0 && y < res_y-1)
							{
								double *kernel_scalar_ptr = &kernel_data[1];
								double *data_scalar_ptr = &i_data[(y-1)*res_x+x];

								data_out += kernel_scalar_ptr[0]*data_scalar_ptr[0];
								data_out += kernel_scalar_ptr[3]*data_scalar_ptr[res_x];
								data_out += kernel_scalar_ptr[6]*data_scalar_ptr[2*res_x];
							}
							else
							{
								for (int j = -1; j <= 1; j++)
								{
									int pos_y = y+j;

									pos_y -= (pos_y >= res_y ? res_y : 0);
									pos_y += (pos_y < 0 ? res_y : 0);

									int idx = (j+1)*3+1;

									double kernel_scalar = kernel_data[idx];
									double data_scalar = i_data[pos_y*res_x+x];

									data_out += kernel_scalar*data_scalar;
								}
							}
						}
					}
					break;

			default:

				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
				for (int y = 0; y < res_y; y++)
				{
					for (int x = 0; x < res_x; x++)
					{
						double &data_out = o_data[y*res_x+x];
						data_out = 0;

						for (int j = -1; j <= 1; j++)
						{
							int pos_y = y+j;

							pos_y -= (pos_y >= res_y ? res_y : 0);
							pos_y += (pos_y < 0 ? res_y : 0);

							assert(pos_y >= 0 && pos_y < res_y);

							for (int i = -1; i <= 1; i++)
							{
								int pos_x = x+i;

								pos_x -= (pos_x >= res_x ? res_x : 0);
								pos_x += (pos_x < 0 ? res_x : 0);

								assert(pos_x >= 0 && pos_x < res_x);

								int idx = (j+1)*3+(i+1);
								assert(idx >= 0 && idx < 9);

								double kernel_scalar = kernel_data[idx];
								double data_scalar = i_data[pos_y*res_x+pos_x];

								data_out += kernel_scalar*data_scalar;
							}
						}
					}
				}
			}
		}
		else
		{
			std::cerr << "Not yet implemented" << std::endl;
		}
	}


	~PlaneData_Kernels()
	{
		if (kernel_data != nullptr)
			MemBlockAlloc::free(kernel_data, sizeof(double)*kernel_size*kernel_size);
	}
};

}

#endif /* SRC_INCLUDE_SWEET_PLANE_PLANEDATA_KERNELS_HPP_ */
