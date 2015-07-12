
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/DataArray.hpp>
#include <sweet/Complex2DArrayFFT.hpp>
#include <sweet/SimulationParameters.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>


SimulationParameters parameters;

int main(int i_argc, char *i_argv[])
{
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);

	SimulationParameters parameters;
	parameters.setup(i_argc, i_argv);

	Complex2DArrayFFT cart(parameters.res);
	Complex2DArrayFFT spec(parameters.res);

	DataArray<2> dataArrayA(parameters.res);

#if 0
	for (std::size_t y = 0; y < dataArrayA.resolution_spec[1]; y++)
		for (std::size_t x = 0; x < dataArrayA.resolution_spec[0]; x++)
//			dataArray.setSpec(y, x, (int)(y+1)*(int)(x+10), -(int)(y+1)*(int)(x+10));
#else
	for (std::size_t y = 0; y < dataArrayA.resolution[1]; y++)
		for (std::size_t x = 0; x < dataArrayA.resolution[0]; x++)
//			dataArrayA.set(y, x, sin(2.0*M_PIl*(double)x/(double)dataArrayA.resolution[0]));
			dataArrayA.set(y, x, sin(2.0*M_PIl*(double)y/(double)dataArrayA.resolution[1]));
#endif

	std::cout << dataArrayA << std::endl;
	std::cout << std::endl;
	dataArrayA.printSpectrum();
	std::cout << std::endl;

	DataArray<2> dataArrayB = dataArrayA.aliasing_scaleUp();
	dataArrayB.setAll(-666);
	dataArrayB = dataArrayA.aliasing_scaleUp();
	dataArrayB.printSpectrum();
	std::cout << std::endl;
	std::cout << dataArrayB << std::endl;


	DataArray<2> dataArrayC = dataArrayB.aliasing_scaleDown(dataArrayA.resolution);
	std::cout << dataArrayC << std::endl;
	dataArrayC.printSpectrum();

	std::cout << "error: " << (dataArrayC-dataArrayA).reduce_sumAbs_quad() << std::endl;

	std::cout << "*************************************************" << std::endl;
	return 1;

	std::cout << "*************************************************" << std::endl;

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
