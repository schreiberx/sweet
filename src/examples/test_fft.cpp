
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/DataArray.hpp>
#include <sweet/SimulationParameters.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>


SimulationParameters parameters;


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


int main(int i_argc, char *i_argv[])
{
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);
	std::cout << std::setprecision(4);
	std::cerr << std::setprecision(4);

	if (i_argc <= 1)
	{
		std::cout << "Please specify spectral or non-spectral test with parameter 1 (spectral) / 0 (non-spectral)" << std::endl;
		return 1;
	}

	SimulationParameters parameters;
	parameters.setup(i_argc, i_argv);


	TestArray cart(parameters.res);
	TestArray spec(parameters.res);

	std::cout << "Cartesian space:" << std::endl;
	cart.setAll(1, 0);
	std::cout << cart << std::endl;

	std::cout << "Spectral space:" << std::endl;
	spec = cart.toSpec();
	std::cout << spec << std::endl;

	std::cout << "*************************************************" << std::endl;

	std::cout << "Cartesian space:" << std::endl;
	cart.setAll(0, 0);
	for (std::size_t y = 0; y < cart.resolution[1]; y++)
		for (std::size_t x = 0; x < cart.resolution[0]; x++)
			cart.set(y, x, pow(-1.0, x+y), 0);
	std::cout << cart << std::endl;

	std::cout << "Spectral space:" << std::endl;
	spec = cart.toSpec();
	std::cout << spec << std::endl;

	std::cout << "*************************************************" << std::endl;

	std::cout << "Cartesian space:" << std::endl;
	cart.setAll(0, 0);
	for (std::size_t y = 0; y < cart.resolution[1]; y++)
		for (std::size_t x = 0; x < cart.resolution[0]; x++)
			cart.set(y, x, std::sin(M_PIl*(x+0.5)/cart.resolution[0]), 0);
	std::cout << cart << std::endl;

	std::cout << "Spectral space:" << std::endl;
	spec = cart.toSpec();
	std::cout << spec << std::endl;

	std::cout << "*************************************************" << std::endl;

	std::cout << "Cartesian space:" << std::endl;
	cart.setAll(0, 0);
	for (std::size_t y = 0; y < cart.resolution[1]; y++)
		for (std::size_t x = 0; x < cart.resolution[0]; x++)
			cart.set(y, x, std::sin(2.0*M_PIl*(x+0.5)/cart.resolution[0]), 0);
	std::cout << cart << std::endl;

	std::cout << "Spectral space:" << std::endl;
	spec = cart.toSpec();
	std::cout << spec << std::endl;

	std::cout << "*************************************************" << std::endl;

	std::cout << "Spectral space:" << std::endl;
	spec.setAll(0, 0);
	for (std::size_t y = 0; y < spec.resolution[1]; y++)
		for (std::size_t x = 0; x < spec.resolution[0]; x++)
			spec.set(y, x, std::sin(2.0*M_PIl*(x+0.5)/spec.resolution[0]), 0);
	std::cout << spec << std::endl;

	std::cout << "Cartesian space:" << std::endl;
	cart = spec.toCart();
	std::cout << cart << std::endl;

	std::cout << "*************************************************" << std::endl;

	std::cout << "Spectral space:" << std::endl;
	spec.setAll(0, 0);
	spec.set(parameters.res[1]-1,parameters.res[0]-1,0,1);
	std::cout << spec << std::endl;

	std::cout << "Cartesian space:" << std::endl;
	cart = spec.toCart();
	std::cout << cart << std::endl;

	std::cout << "*************************************************" << std::endl;

	std::cout << "Cartesian space:" << std::endl;
	cart.setAll(0, 0);
	for (std::size_t y = 0; y < cart.resolution[1]; y++)
		for (std::size_t x = 0; x < cart.resolution[0]; x++)
			cart.set(y, x, (2.0*M_PIl*(x)/cart.resolution[0]), 0);
	std::cout << cart << std::endl;

	std::cout << "Spectral space:" << std::endl;
	spec = cart.toSpec();
	std::cout << spec << std::endl;

	std::cout << "*************************************************" << std::endl;

	fftw_cleanup();
	exit(-1);

	return 1;
}
