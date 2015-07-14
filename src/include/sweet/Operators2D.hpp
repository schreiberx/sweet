/*
 * Operators.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_OPERATORS2D_HPP_
#define SRC_INCLUDE_SWEET_OPERATORS2D_HPP_


#if SWEET_USE_SPECTRAL_SPACE
	#include "Complex2DArrayFFT.hpp"
#endif



class Operators2D
{
public:
	// differential operators (central / forward / backward)
	DataArray<2> diff_c_x, diff_c_y;
	DataArray<2> diff_f_x, diff_f_y;
	DataArray<2> diff_b_x, diff_b_y;

	DataArray<2> diff2_c_x, diff2_c_y;

	DataArray<2> avg_f_x, avg_f_y;
	DataArray<2> avg_b_x, avg_b_y;

	DataArray<2> shift_left;
	DataArray<2> shift_right;
	DataArray<2> shift_up;
	DataArray<2> shift_down;

	inline DataArray<2> arakawa_jacobian(
			const DataArray<2> &i_a,
			const DataArray<2> &i_b
	)
	{
//		return diff_c_x(i_a)*diff_c_y(i_b) - diff_c_y(i_a)*diff_c_x(i_b);
		return diff_c_y(i_a) - diff_c_x(i_a);
	}



	/**
	 *        __2
	 * apply  \/  operator (aka Laplace)
	 */
	inline DataArray<2> laplace(
			const DataArray<2> &i_a
	)
	{
		return diff2_c_x(i_a)+diff2_c_y(i_a);
	}


	/**
	 *        __
	 * apply  \/ .  operator
	 */
	inline DataArray<2> diff_dot(
			const DataArray<2> &i_a
	)
	{
		return diff_c_x(i_a)+diff_c_y(i_a);
	}

	Operators2D(
		std::size_t res[2],		///< resolution
		double i_domain_size[2],	///< domain size
		bool i_use_spectral_diffs = false
	)	:
		diff_c_x(res),
		diff_c_y(res),

		diff_f_x(res),
		diff_f_y(res),
		diff_b_x(res),
		diff_b_y(res),

		diff2_c_x(res),
		diff2_c_y(res),

		avg_f_x(res),
		avg_f_y(res),
		avg_b_x(res),
		avg_b_y(res),

		shift_left(res),
		shift_right(res),
		shift_up(res),
		shift_down(res)
	{
		double h[2] = {(double)i_domain_size[0] / (double)res[0], (double)i_domain_size[1] / (double)res[1]};

		/////////////////////////////////////////////////////////////////////

		double avg_f_x_kernel[3][3] = {
				{0,0,0},
				{0,1,1},
				{0,0,0},
		};
		avg_f_x.setup_kernel(avg_f_x_kernel, 0.5);

		double avg_f_y_kernel[3][3] = {
				{0,1,0},
				{0,1,0},
				{0,0,0},
		};
		avg_f_y.setup_kernel(avg_f_y_kernel, 0.5);

		double avg_b_x_kernel[3][3] = {
				{0,0,0},
				{1,1,0},
				{0,0,0},
		};
		avg_b_x.setup_kernel(avg_b_x_kernel, 0.5);

		double avg_b_y_kernel[3][3] = {
				{0,0,0},
				{0,1,0},
				{0,1,0},
		};
		avg_b_y.setup_kernel(avg_b_y_kernel, 0.5);



		/////////////////////////////////////////////////////////////////////

		if (i_use_spectral_diffs)
		{
#if !SWEET_USE_SPECTRAL_SPACE
			std::cerr << "Activate FFTW during compile time to use spectral diffs" << std::endl;
			exit(-1);
#else

			std::size_t *res = diff2_c_x.resolution;

			Complex2DArrayFFT spec_dx(res);

			for (int j = 0; j < (int)res[1]; j++)
			{
				for (int i = 0; i < (int)res[0]; i++)
				{
					if (i < (int)res[0]/2)
						spec_dx.set(j, i, 0, (double)i*(2.0*M_PIl)/(double)i_domain_size[0]);
					else if (i == (int)res[0]/2)
						spec_dx.set(j, i, 0, 0);
					else
						spec_dx.set(j, i, 0, (-(double)res[0]+(double)i)*(2.0*M_PIl)/(double)i_domain_size[0]);
				}
			}

			Complex2DArrayFFT cart_x = spec_dx.toCart();

			double normalize = 1.0/(diff_c_x.resolution[0]*diff_c_x.resolution[1]);

			for (int j = 0; j < (int)diff_c_x.resolution[1]; j++)
					for (int i = 0; i < (int)diff_c_x.resolution[0]; i++)
						diff_c_x.set(j, i, cart_x.getRe(j, i)*normalize);

			Complex2DArrayFFT spec_dy(res);

			for (int i = 0; i < (int)res[0]; i++)
			{
				for (int j = 0; j < (int)res[1]; j++)
				{
					if (j < (int)res[1]/2)
						spec_dy.set(j, i, 0, (double)j*(2.0*M_PIl)/(double)i_domain_size[1]);
					else if (j == (int)res[1]/2)
						spec_dy.set(j, i, 0, 0);
					else
						spec_dy.set(j, i, 0, (-(double)res[1]+(double)j)*(2.0*M_PIl)/(double)i_domain_size[1]);
				}
			}

			Complex2DArrayFFT cart_y = spec_dy.toCart();

			for (int j = 0; j < (int)diff_c_y.resolution[1]; j++)
					for (int i = 0; i < (int)diff_c_y.resolution[0]; i++)
						diff_c_y.set(j, i, cart_y.getRe(j, i)*normalize);


			diff2_c_x = diff_c_x(diff_c_x);
			diff2_c_y = diff_c_y(diff_c_y);

#endif
		}
		else
		{
			double diff1_x_kernel[3][3] = {
					{0,0,0},
					{-1.0,0,1.0},
					{0,0,0}
			};
			diff_c_x.setup_kernel(diff1_x_kernel, 1.0/(2.0*h[0]));

			double diff1_y_kernel[3][3] = {
					{0,1.0,0},	// higher y coordinate
					{0,0,0},
					{0,-1.0,0},	// lower y coordinate
			};
			diff_c_y.setup_kernel(diff1_y_kernel, 1.0/(2.0*h[1]));


			double d_f_x_kernel[3][3] = {
					{0,0,0},
					{0,-1,1},
					{0,0,0}
			};
			diff_f_x.setup_kernel(d_f_x_kernel, 1.0/h[0]);

			double d_f_y_kernel[3][3] = {
					{0,1,0},
					{0,-1,0},
					{0,0,0},
			};
			diff_f_y.setup_kernel(d_f_y_kernel, 1.0/h[1]);


			double d_b_x_kernel[3][3] = {
					{0,0,0},
					{-1,1,0},
					{0,0,0}
			};
			diff_b_x.setup_kernel(d_b_x_kernel, 1.0/h[0]);

			double d_b_y_kernel[3][3] = {
					{0,0,0},
					{0,1,0},
					{0,-1,0},
			};
			diff_b_y.setup_kernel(d_b_y_kernel, 1.0/h[1]);


			double diff2_x_kernel[3][3] = {
					{0,0,0},
					{1.0,-2.0,1.0},
					{0,0,0}
				};
			diff2_c_x.setup_kernel(diff2_x_kernel, 1.0/(h[0]*h[0]));

			double diff2_y_kernel[3][3] = {
					{0,1.0,0},
					{0,-2.0,0},
					{0,1.0,0}
			};
			diff2_c_y.setup_kernel(diff2_y_kernel, 1.0/(h[1]*h[1]));
		}

		/////////////////////////////////////////////////////////////////////

		double shift_left_kernel[3][3] = {
				{0,0,0},
				{0,0,1},
				{0,0,0},
		};
		shift_left.setup_kernel(shift_left_kernel);

		double shift_right_kernel[3][3] = {
				{0,0,0},
				{1,0,0},
				{0,0,0},
		};
		shift_right.setup_kernel(shift_right_kernel);

		double shift_up_kernel[3][3] = {
				{0,0,0},
				{0,0,0},
				{0,1,0},
		};
		shift_up.setup_kernel(shift_up_kernel);

		double shift_down_kernel[3][3] = {
				{0,1,0},
				{0,0,0},
				{0,0,0},
		};
		shift_down.setup_kernel(shift_down_kernel);

	}
};



#endif /* SRC_INCLUDE_SWEET_OPERATORS2D_HPP_ */
