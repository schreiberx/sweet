
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/DataArray.hpp>
#include <sweet/SimulationParameters.hpp>
#include <sweet/Operators2D.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>

SimulationParameters parameters;

#if 0
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


int main(int i_argc, char *i_argv[])
{
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);
	std::cout << std::setprecision(4);
	std::cerr << std::setprecision(4);

	SimulationParameters parameters;
	parameters.use_spectral_diffs = 1;
	parameters.setup(i_argc, i_argv);

	if (parameters.use_spectral_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;

	double prev_error_x = 0;
	double prev_error_y = 0;

	double freq_x = 4.0;
	double freq_y = 4.0;

#if 0
	{
		TestArray spec(parameters.res);

		spec.setAll(123, 456);
		for (int j = 0; j < parameters.res[1]; j++)
			for (int i = 0; i < parameters.res[0]; i++)
			{
				if (i < parameters.res[0]/2)
					spec.set(j, i, 0, i);
				else if (i == parameters.res[0]/2)
					spec.set(j, i, 0, 0);
				else
					spec.set(j, i, 0, -(int)parameters.res[0]+(int)i);
			}
		std::cout << spec << std::endl;
		std::cout << std::endl;

		TestArray cart = spec.toCart();
		std::cout << cart << std::endl;
		exit(1);
	}
#endif

	for (std::size_t res_x = parameters.res[0]; res_x <= 4096; res_x *= 2)
//	for (std::size_t res_y = 32; res_y <= 4096; res_y *= 2)
	{
#if 0
		if (res_x != parameters.res[0])
			break;
#endif
		std::size_t res_y = res_x;
		std::size_t res[2] = {res_x,res_y};

		parameters.res[0] = res[0];
		parameters.res[1] = res[1];
		parameters.reset();

		DataArray<2> h(res);
		DataArray<2> h_diff_x(res);
		DataArray<2> h_diff_y(res);

		Operators2D op(parameters.res, parameters.sim_domain_length, parameters.use_spectral_diffs);

		// scale factor to avoid very very tiny values
		double scale = ((double)parameters.sim_domain_length[0]*(double)parameters.sim_domain_length[1]);
//		scale *= scale;
//		double scale = 1.0;

		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				double x = ((double)i)/(double)parameters.res[0];
				double y = ((double)j)/(double)parameters.res[1];

				h.set(
					j, i,
					sin(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)*scale
				);

				h_diff_x.set(
					j, i,
					freq_x*M_PIl*cos(freq_x*M_PIl*x)*cos(freq_y*M_PIl*y)/(double)parameters.sim_domain_length[0]*scale
				);

				h_diff_y.set(
					j, i,
					sin(freq_x*M_PIl*x)*freq_y*M_PIl*(-sin(freq_y*M_PIl*y))/(double)parameters.sim_domain_length[1]*scale
				);
			}
		}

		double normalization = 1.0/((double)res[0]*(double)res[1]);
		double err_x = (op.diff_c_x(h)-h_diff_x).reduce_sumAbs()*normalization;
		double err_y = (op.diff_c_y(h)-h_diff_y).reduce_sumAbs()*normalization;

		if (parameters.use_spectral_diffs)
		{
			std::cout << res_x << "\t" << res_y << "\t" << err_x << "\t" << err_y << std::endl;

			if (err_x > 10e-10)
				std::cerr << "Error threshold for diff-X too high for spectral differentiation!" << std::endl;

			if (err_y > 10e-10)
				std::cerr << "Error threshold for diff-Y too high for spectral differentiation!" << std::endl;
		}
		else
		{
			std::cout << res_x << "\t" << res_y << "\t" << err_x << "\t" << err_y << "\t" << prev_error_x/err_x << "\t" << prev_error_y/err_y << std::endl;

			prev_error_x = err_x;
			prev_error_y = err_y;
		}
	}


	return 1;
}
