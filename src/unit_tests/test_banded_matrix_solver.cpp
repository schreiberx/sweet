/*
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#include <iostream>
#include <iomanip>

#include <libmath/BandedMatrixPhysicalReal.hpp>
#include <libmath/BandedMatrixSolver.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>

#include <sweet/SWEETError.hpp>

#define T std::complex<double>

void getRandomRHS(
			SphereData_Spectral& o_rhs
		)
{
	int max_random = 1;
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<> distr(0, max_random); // define the range
	for (int i = 0; i < o_rhs.sphereDataConfig->spectral_array_data_number_of_elements; i++)
		o_rhs.spectral_space_data[i] = (T)distr(gen) / (T)max_random;
}


void getRandomSystem(
			int i_size,
			int i_ndiag,
			BandedMatrixSolver< std::complex<double> >& i_bandedMatrixSolver,
			T* o_lhs,
			T* o_rhs
		)
{
	int max_random = 1000;
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<> distr(0, max_random); // define the range
	for (int i = 0; i < i_size; i++)
	{
		o_rhs[i] = (T)distr(gen) / (T)max_random;
		for (int j = std::max(0, i - i_ndiag); j < std::min(i_size, i + i_ndiag + 1); j++)
		{
			int idx = i_bandedMatrixSolver.getMatrixIndex(i_size, i, j);
			o_lhs[idx] = (T)distr(gen) / (T)max_random;
			std::cout << "SETTING " << i << " " << j << " " << idx << " " << i_size * i_bandedMatrixSolver.num_diagonals << " " << o_lhs[idx] << std::endl;
		}
	}
}


/////////int main(
/////////		int i_argc,
/////////		char *const i_argv[]
/////////)
/////////{
/////////
/////////	// Test banded matrix solver for given matrix sizes (N = 2^n) and numbers of outer diagonals (ndiag U + ndiag L)
/////////	for (int n = 0; n < 12; n++)
/////////	{
/////////		int N = int(std::pow(2, n));
/////////		SphereData_Config* sphereDataConfig;
/////////
/////////		for (int ndiag = 0; ndiag < 10; ndiag++)
/////////		{
/////////
/////////			int size = sphereDataConfig->spectral_modes_n_max + 1;
/////////
/////////			BandedMatrixPhysicalReal<T> lhs;
/////////			BandedMatrixSolver< std::complex<double> > bandedMatrixSolver;
/////////
/////////			SphereData_Spectral rhs(sphereDataConfig);
/////////			SphereData_Spectral x(sphereDataConfig);
/////////
/////////			lhs.setup(sphereDataConfig, ndiag);
/////////			lhs.setRandomElements();
/////////			getRandomRHS(rhs);
/////////			bandedMatrixSolver.setup(size, ndiag);
/////////
/////////
/////////			for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
/////////			{
/////////
/////////				double time_pivoting;
/////////				double time_solving;
/////////
/////////				int idx = sphereDataConfig->getArrayIndexByModes(m,m);
/////////
/////////				bandedMatrixSolver.solve_diagBandedInverse(
/////////									sphereDataConfig->spectral_modes_n_max+1-m,					// size of block
/////////									&lhs.data[idx*lhs.num_diagonals],						// A
/////////									&rhs.spectral_space_data[idx],							// b
/////////									&x.spectral_space_data[idx],							// x
/////////									time_pivoting,
/////////									time_solving
/////////									);
/////////
/////////				bool ok = bandedMatrixSolver.checkSolution(
/////////									sphereDataConfig->spectral_modes_n_max+1-m,					// size of block
/////////									&lhs.data[idx*lhs.num_diagonals],						// A
/////////									&rhs.spectral_space_data[idx],							// b
/////////									&x.spectral_space_data[idx]							// x
/////////								);
/////////
/////////				if (!ok)
/////////					SWEETError("Wrong system solution for N = " + std::to_string(n) + " and ndiag = " + std::to_string(ndiag));
/////////			}
/////////
/////////		}
/////////	}
/////////
/////////}


int main(
		int i_argc,
		char *const i_argv[]
)
{

	// Test banded matrix solver for given matrix sizes (N = 2^n) and numbers of outer diagonals (ndiag U + ndiag L)
	for (int n = 0; n < 12; n++)
	{
		int size = int(std::pow(2, n));

		for (int ndiag = 0; ndiag < 10; ndiag++)
		{

			if (ndiag + 1 > size)
				continue;

			std::cout << std::endl << std::endl;
			std::cout << "SOLVING FOR size,ndiag = " << size << "," << ndiag << std::endl;

			BandedMatrixSolver< std::complex<double> > bandedMatrixSolver;
			bandedMatrixSolver.setup(size, ndiag);

			T* lhs = nullptr;
			T* rhs = nullptr;
			T* x = nullptr;
			lhs = MemBlockAlloc::alloc<T>(size * (2 * ndiag + 1) * sizeof(T));
			rhs = MemBlockAlloc::alloc<T>(size * sizeof(T));
			x = MemBlockAlloc::alloc<T>(size * sizeof(T));

			getRandomSystem(size, ndiag, bandedMatrixSolver, lhs, rhs);

			double time_pivoting = 0;
			double time_solving = 0;

			bandedMatrixSolver.printAugmentedMatrix(size, lhs, rhs);

			bandedMatrixSolver.solve_diagBandedInverse(
								size,					// size of block
								lhs,						// A
								rhs,							// b
								x,							// x
								time_pivoting,
								time_solving
								);

			bool ok = bandedMatrixSolver.checkSolution(
								size,					// size of block
								lhs,						// A
								rhs,							// b
								x							// x
							);

			if (!ok)
				SWEETError("Wrong system solution for N = " + std::to_string(n) + " and ndiag = " + std::to_string(ndiag));

			MemBlockAlloc::free(lhs, size * ndiag * sizeof(T));
			MemBlockAlloc::free(rhs, size * sizeof(T));
			MemBlockAlloc::free(x, size * sizeof(T));
		}
	}

}
