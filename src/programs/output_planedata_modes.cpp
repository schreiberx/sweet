/*
 * output_plane_data_modes.cpp
 *
 *  Created on: 18 Juli 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */


#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperators.hpp>

SimulationVariables simVars;

PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

class PlaneDataModes
{
public:
	PlaneOperators op;


	void output_planedata()
	{
		PlaneData test(planeDataConfig);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				test.physical_set_zero();
				test.p_physical_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Phys coord : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalArrayData();
				std::cout << "***" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(1.0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}
	}



	void output_planedatacomplex()
	{
		PlaneDataComplex test(planeDataConfig);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				test.physical_set_zero();
				test.p_physical_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Phys coord : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_complex_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_complex_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_complex_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_complex_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_complex_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_complex_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(1.0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}
	}

	PlaneDataModes()	:
		op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
	{
		output_planedata();
		output_planedatacomplex();
	}
};



int main(
		int i_argc,
		char *i_argv[]
)
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	simVars.outputConfig();

	PlaneDataModes planeDataModes;

	return 0;
}
