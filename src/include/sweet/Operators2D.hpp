/*
 * Operators.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_OPERATORS2D_HPP_
#define SRC_INCLUDE_SWEET_OPERATORS2D_HPP_


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

	Operators2D(
		double h[2],			///< cell sizes,
		std::size_t res[2]		///< resolution
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
