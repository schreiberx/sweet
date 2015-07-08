/*
 * Operators.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_OPERATORS2D_HPP_
#define SRC_INCLUDE_SWEET_OPERATORS2D_HPP_


#if SWEET_USE_SPECTRAL_SPACE
class TestArray
{
public:
	std::size_t resolution[2];

	fftw_plan plan_to_cart;
	fftw_plan plan_to_spec;

	double *data;

	void setup_fftw()
	{
		plan_to_spec =
				fftw_plan_dft_2d(
					resolution[1],	// n0 = ny
					resolution[0],	// n1 = nx
					(fftw_complex*)data,
					(fftw_complex*)data,
					FFTW_FORWARD,
					0
				);


		plan_to_cart =
				fftw_plan_dft_2d(
					resolution[1],	// n0 = ny
					resolution[0],	// n1 = nx
					(fftw_complex*)data,
					(fftw_complex*)data,
					FFTW_BACKWARD,
					0
				);

		if (plan_to_spec == nullptr)
		{
			std::cerr << "Failed to create plan_backward for fftw" << std::endl;
			exit(-1);
		}
	}

public:
	TestArray(
			std::size_t i_res[2]
	)	:
		plan_to_cart(nullptr),
		plan_to_spec(nullptr)
	{
		resolution[0] = i_res[0];
		resolution[1] = i_res[1];

		data = alloc_aligned_mem<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		setup_fftw();
	}


public:
	TestArray(const TestArray &i_testArray)
	:
		plan_to_cart(nullptr),
		plan_to_spec(nullptr)
	{
		resolution[0] = i_testArray.resolution[0];
		resolution[1] = i_testArray.resolution[1];

		data = alloc_aligned_mem<double>(sizeof(double)*resolution[0]*resolution[1]*2);

		setup_fftw();

		memcpy(data, i_testArray.data, resolution[0]*resolution[1]*2);
	}


public:
	TestArray& operator=(const TestArray &i_testArray)
	{
		resolution[0] = i_testArray.resolution[0];
		resolution[1] = i_testArray.resolution[1];

		memcpy(data, i_testArray.data, resolution[0]*resolution[1]*2*sizeof(double));
		return *this;
	}

	~TestArray()
	{
		fftw_free(plan_to_spec);

		free(data);
	}


	TestArray toSpec()
	{
		TestArray o_testArray(resolution);

		fftw_execute_dft(
				plan_to_spec,
				(fftw_complex*)this->data,
				(fftw_complex*)o_testArray.data
			);


		return o_testArray;
	}


	TestArray toCart()
	{
		TestArray o_testArray(resolution);

		fftw_execute_dft(
				plan_to_cart,
				(fftw_complex*)this->data,
				(fftw_complex*)o_testArray.data
			);

		return o_testArray;
	}


	void set(int y, int x, double re, double im)
	{
		data[(y*resolution[0]+x)*2+0] = re;
		data[(y*resolution[0]+x)*2+1] = im;
	}

	double getRe(int y, int x)	const
	{
		return data[(y*resolution[0]+x)*2+0];
	}

	double getIm(int y, int x)	const
	{
		return data[(y*resolution[0]+x)*2+1];
	}

	void setAll(double re, double im)
	{
		for (std::size_t y = 0; y < resolution[1]; y++)
			for (std::size_t x = 0; x < resolution[0]; x++)
				set(y, x, re, im);
	}


	friend
	inline
	std::ostream& operator<<(std::ostream &o_ostream, const TestArray &i_testArray)
	{
		for (int y = i_testArray.resolution[1]-1; y >= 0; y--)
		{
			for (std::size_t x = 0; x < i_testArray.resolution[0]; x++)
			{
				double value_re = i_testArray.getRe(y, x);
				double value_im = i_testArray.getIm(y, x);
				o_ostream << "(" << value_re << ", " << value_im << ")\t";
			}
			o_ostream << std::endl;
		}
		return o_ostream;
	}
};
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
			double wx = (2.0*M_PIl)/(double)i_domain_size[0];
			double wy = (2.0*M_PIl)/(double)i_domain_size[1];
#if 0
			// See
			// http://stackoverflow.com/questions/16749349/fftw-simple-derivative-in-fortran

			// DIFF X
			for (int j = 0; j < (int)diff_c_x.resolution_spec[1]; j++)
				for (int i = 0; i < (int)diff_c_x.resolution_spec[0]-1; i++)
					diff_c_x.setSpec(j, i, 0, (double)i*wx);

			for (int j = 0; j < (int)diff_c_x.resolution_spec[1]; j++)
				diff_c_x.setSpec(j, (int)diff_c_x.resolution_spec[0]-1, 0, 0);

#else

			TestArray spec_dx(res);

	//				spec_dx.setAll(123, 456);
			for (int j = 0; j < (int)res[1]; j++)
				for (int i = 0; i < (int)res[0]; i++)
				{
					if (i < (int)res[0]/2)
						spec_dx.set(j, i, 0, i*wx);
					else if (i == (int)res[0]/2)
						spec_dx.set(j, i, 0, 0);
					else
						spec_dx.set(j, i, 0, (-(int)res[0]+(int)i)*wx);
				}

				TestArray cart_x = spec_dx.toCart();

				double normalize = 1.0/(diff_c_y.resolution[0]*diff_c_y.resolution[1]);
				for (int j = 0; j < (int)diff_c_y.resolution[1]; j++)
						for (int i = 0; i < (int)diff_c_y.resolution[0]; i++)
							diff_c_x.set(j, i, cart_x.getRe(j, i)*normalize);

#endif

#if 0
			// DIFF Y
			for (int j = 0; j < ((int)diff_c_y.resolution_spec[1])/2-1; j++)
			{
				for (int i = 0; i < (int)diff_c_y.resolution_spec[0]; i++)
					diff_c_y.setSpec(j, i, 0, (double)j*wy);

				// last row = 0
				diff_c_y.setSpec(j, diff_c_y.resolution_spec[0]-1, 0,0);
			}

			for (int j = (int)diff_c_y.resolution_spec[1]/2; j < (int)diff_c_y.resolution_spec[1]; j++)
			{
				for (int i = 0; i < (int)diff_c_y.resolution_spec[0]-1; i++)
					diff_c_y.setSpec(j, i, 0, ((double)j-(double)diff_c_y.resolution_spec[1])*wy);

				// last row = 0
				diff_c_y.setSpec(j, diff_c_y.resolution_spec[0]-1, 0,0);
			}
#else
			TestArray spec_dy(res);

//				spec_dy.setAll(123, 456);
			for (int i = 0; i < (int)res[0]; i++)
				for (int j = 0; j < (int)res[1]; j++)
				{
					if (j < (int)res[1]/2)
						spec_dy.set(j, i, 0, j*wy);
					else if (j == (int)res[1]/2)
						spec_dy.set(j, i, 0, 0);
					else
						spec_dy.set(j, i, 0, (-(int)res[1]+(int)j)*wy);
				}

				TestArray cart_y = spec_dy.toCart();

				for (int j = 0; j < (int)diff_c_y.resolution[1]; j++)
						for (int i = 0; i < (int)diff_c_y.resolution[0]; i++)
							diff_c_y.set(j, i, cart_y.getRe(j, i)*normalize);
#endif

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
