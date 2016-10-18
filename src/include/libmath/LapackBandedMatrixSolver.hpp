/*
 * DiagBandedMatrix.hpp
 *
 *  Created on: 24 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_LIBMATH_LAPACKBANDEDMATRIXSOLVER_HPP_
#define SRC_INCLUDE_LIBMATH_LAPACKBANDEDMATRIXSOLVER_HPP_

#include <complex>
#include <string.h>
#include <stdlib.h>


/*
 * LAPACK description:
 *
 * For an array
 * 	REAL A( LDA, * )
 *
 * LDA is called the leading dimension of the array.
 * In Fortran, the values in this dimension
 * are consecutively stored in memory.
 *
 * aij = A(i,j)
 *
 *    a11 a12 a13 a14   |
 *    a21 a22 a23 a24   |
 *    a31 a32 a33 a34   |LDA
 *    a41 a42 a43 a44   v
 *
 */
/*
 * We can't store the inverse of the general band matrix since
 * the inverse would result in a dense matrix.
 *
 * http://www.netlib.org/lapack/explore-html/d9/dbb/group__complex16_g_bsolve_ga908abc0aad64131b9a32edb08510eb00.html#ga908abc0aad64131b9a32edb08510eb00
 *
 * ZGBSV computes the solution to system of linear equations A * X = B for GB matrices (simple driver)
 *
subroutine zgbsv
	(
		integer  	N,
		integer  	KL,
		integer  	KU,
		integer  	NRHS,
		complex*16, dimension( ldab, * )  	AB,
		integer  	LDAB,
		integer, dimension( * )  	IPIV,
		complex*16, dimension( ldb, * )  	B,
		integer  	LDB,
		integer  	INFO
	)
*/

/*
 * We create the interface to this function here
 */
extern "C"
{
	void zgbsv_(
			const int &N,
			const int &KL,
			const int &KU,
			const int &NRHS,
			std::complex<double> *AB,
			const int &LDAB,
			int *IPIV,
			std::complex<double> *B,
			const int &LDB,
			int &INFO
	);
#if 0
	void zlapmr_(
			int &forward,
			int &M,	// rows
			int &N,	// cols
			std::complex<double> *X,
			int &LDX,	// leading dimension
			int *K		// pivotization table
	);
#endif
}


template <typename T>
class BandedMatrixSolverCommon
{
public:
	int max_N;
	int num_diagonals;
	int num_halo_size_diagonals;

	int LDAB;

	std::complex<double>* AB;
	int *IPIV;

	BandedMatrixSolverCommon()	:
		AB(nullptr),
		IPIV(nullptr)
	{
	}



	~BandedMatrixSolverCommon()
	{
		shutdown();
	}

	void setup(
			int i_max_N,			///< size of the matrix
			int i_num_off_diagonals		///< number of block diagonals
	)
	{
		max_N = i_max_N;
		num_diagonals = 2*i_num_off_diagonals+1;
		num_halo_size_diagonals = i_num_off_diagonals;

		assert(2*num_halo_size_diagonals+1 == num_diagonals);

		LDAB = 2*num_halo_size_diagonals + num_halo_size_diagonals + 1;

		AB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*LDAB*i_max_N);
		IPIV = (int*)malloc(sizeof(int)*i_max_N);
	}


	void shutdown()
	{
		if (IPIV != nullptr)
		{
			free(IPIV);
			IPIV = nullptr;
		}

		if (AB != nullptr)
		{
			free(AB);
			AB = nullptr;
		}
	}


	void print_array_fortran(
			const T *i_data,
			int i_cols,
			int i_rows
	)
	{
		// rows
		for (int j = 0; j < i_rows; j++)
		{
			std::cout << j << ": ";
			// cols
			for (int i = 0; i < i_cols; i++)
			{
				std::cout << i_data[j+i*i_rows] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}


	void print_array_c(
			const T *i_data,
			int i_cols,
			int i_rows
	)
	{
		// rows
		for (int j = 0; j < i_rows; j++)
		{
			std::cout << j << ": ";
			// cols
			for (int i = 0; i < i_cols; i++)
			{
				std::cout << i_data[j*i_cols+i] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
};

template <typename T>
class SphBandedMatrix	:
		public BandedMatrixSolverCommon<T>
{
public:
	/**
	 * Solve a diagonal banded matrix
	 *
	 * A*X = B
	 */
	void solve_diagBandedInverse_C(
		const std::complex<double>* i_A,		///< Matrix for input and in-place transformations
		const std::complex<double>* i_b,		///< RHS of equation and output of solution X
		std::complex<double>* o_x,
		int i_size
	);


	void solve_diagBandedInverse_Fortran(
		const std::complex<double>* i_A,		///< Matrix for input and in-place transformations
		const std::complex<double>* i_b,		///< RHS of equation and output of solution X
		std::complex<double>* o_x,
		int i_size
	);
};


template <>
class SphBandedMatrix<std::complex<double>>	:
		public BandedMatrixSolverCommon<std::complex<double>>
{
	typedef std::complex<double> T;

	/**
	 * Solve for input matrix
	 *
	 * i_A: cols: num_diagonals
	 *      rows: i_size
	 *
	 * The Fortran array will be transposed and with a size of (rows: LDAB, cols: i_size)
	 *
	 * i_b: RHS of equation
	 *
	 * o_x: Solution
	 */
public:
	void solve_diagBandedInverse_Carray(
		const std::complex<double>* i_A,
		const std::complex<double>* i_b,
		std::complex<double>* o_x,
		int i_size
	)
	{
		assert(max_N >= i_size);

#if 0
		/**
		 * WARNING: LEAVE THIS BLOCK FOR DEBUGGING PURPOSE!!!
		 *
		 * 1) Decodes the compact stored matrix A into a full NxN matrix.
		 * 2) Convert to Fortran storage format
		 * 3) convert to LAPACK general band matrix format
		 */
		std::cout << "C ARRAY: A compactified" << std::endl;
		print_array_c(i_A, num_diagonals, i_size);

		T *fortran_A = new T[i_size*i_size];

		for (int i = 0; i < i_size*i_size; i++)
			fortran_A[i] = std::numeric_limits<double>::infinity();

		// source c cols
		for (int i = 0; i < num_diagonals; i++)
		{
			// source c rows
			for (int j = 0; j < i_size; j++)
			{
				int di = j+i-num_halo_size_diagonals;

				if (di < 0 || di >= i_size)
					continue;

				int dj = j;
				fortran_A[di*i_size+dj] = i_A[j*num_diagonals+i];
			}
		}

		std::cout << "FORTRAN ARRAY: A expanded" << std::endl;
		print_array_fortran(fortran_A, i_size, i_size);

		for (int i = 0; i < i_size*LDAB; i++)
			AB[i] = std::numeric_limits<double>::infinity();

		for (int i = 0; i < i_size; i++)
		{
			int fi = i+1;

			// columns for output fortran array
			for (int j = 0; j < i_size; j++)
			{
				int fj = j+1;

				if (std::max(1, fj-num_halo_size_diagonals) <= fi && fi <= std::min(i_size, fj+num_halo_size_diagonals))
					AB[(num_diagonals+fi-fj)-1 + (fj-1)*LDAB] = fortran_A[i+j*i_size];
			}
		}


		std::cout << "FORTRAN ARRAY: A compactified" << std::endl;
		print_array_fortran(AB, i_size, LDAB);


#if 1
		for (int i = 0; i < i_size*LDAB; i++)
			AB[i] = std::numeric_limits<double>::infinity();

		// columns for output fortran array
		// rows for input c array
		for (int j = 0; j < i_size; j++)
		{
			// rows for output fortran array
			// columns for input c array
			for (int i = 0; i < num_diagonals; i++)
			{
				// compute square matrix indices
				int si = j+(num_halo_size_diagonals-i);
				int sj = j;

//				std::cout << sj << " " << si << std::endl;

				if (si < 0 || si >= i_size)
					continue;


				// AB is LDAB large!
				assert(LDAB*max_N > i*i_size+j);
				assert(LDAB*max_N > i+j*num_diagonals);

				AB[(num_diagonals+si-sj-1) + sj*LDAB] = i_A[(j-i+num_halo_size_diagonals)*num_diagonals + i];
			}
		}
#endif

		delete fortran_A;

		std::cout << "FORTRAN ARRAY: A compactified (alternative, should match previous one)" << std::endl;
		print_array_fortran(AB, i_size, LDAB);

#else

#ifndef NDEBUG
		for (int i = 0; i < i_size*LDAB; i++)
			AB[i] = std::numeric_limits<double>::infinity();
#endif

		// columns for output fortran array
		// rows for input c array
		for (int j = 0; j < i_size; j++)
		{
			// rows for output fortran array
			// columns for input c array
			for (int i = 0; i < num_diagonals; i++)
			{
				// compute square matrix indices
				int si = j+(num_halo_size_diagonals-i);
				int sj = j;

//				std::cout << sj << " " << si << std::endl;

				if (si < 0 || si >= i_size)
					continue;


				// AB is LDAB large!
				assert(LDAB*max_N > i*i_size+j);
				assert(LDAB*max_N > i+j*num_diagonals);

				AB[(num_diagonals+si-sj-1) + sj*LDAB] = i_A[(j-i+num_halo_size_diagonals)*num_diagonals + i];
			}
		}
#endif

		solve_diagBandedInverse_FortranArray(AB, i_b, o_x, i_size);
	}



public:
	void solve_diagBandedInverse_FortranArray(
		const std::complex<double>* i_A,	///< A of max size
		const std::complex<double>* i_b,
		std::complex<double>* o_x,
		int i_size
	)
	{
		/*
		 * Make a copy of the array data since this is a destructive function
		 */
		if (AB != i_A)
			memcpy((void*)AB, (const void*)i_A, sizeof(std::complex<double>)*num_diagonals*LDAB);

		memcpy((void*)o_x, (const void*)i_b, sizeof(std::complex<double>)*i_size);

		solve_diagBandedInverse_FortranArray_inplace(AB, o_x, i_size);
	}



public:
	void solve_diagBandedInverse_FortranArray_inplace(
		std::complex<double>* io_A,		///< A of max size
		std::complex<double>* io_b_x,	///< rhs and solution x
		int i_size
	)
	{
		assert(num_diagonals & 1 == 1);
		assert(AB != nullptr);

		int info;

#if 0
		std::cout << "************************************" << std::endl;
		std::cout << "i_size: " << i_size << std::endl;
		std::cout << "num_halo_size_diagonals: " << num_halo_size_diagonals << std::endl;
		std::cout << "LDAB: " << LDAB << std::endl;
#endif

		zgbsv_(
				i_size,				// number of linear equations
				num_halo_size_diagonals,	// number of subdiagonals
				num_halo_size_diagonals,	// number of superdiagonals
				1,				// number of columns of matrix B
				io_A,				// array with matrix A to solve for
				LDAB,				// leading dimension of matrix A
				IPIV,				// integer array for pivoting
				io_b_x,				// output array
				i_size,				// leading dimension of array o_x
				info
			);

		if (info != 0)
		{
			std::cerr << "zgbsv returned INFO != 0: " << info << std::endl;
			assert(false);
			exit(1);
		}

#if 0
		std::cout << std::endl;
		for (int i = 0; i < i_size; i++)
			std::cout << IPIV[i] << ", ";
		std::cout << std::endl;

		std::cout << std::endl;
		for (int i = 0; i < i_size; i++)
			std::cout << o_x[i] << ", ";
		std::cout << std::endl;
#endif

#if 0
		/**
		 * TODO: This shouldn't be required since the solution is directly computed.
		 *
		 * TODO: Check pivotization
		 */
		int bvalue = true;
		zlapmr_(
				bvalue,	// true = forward permutation
				one,	// rows
				i_size,	// cols
				o_x,	// data
				i_size,	// leading dimension
				IPIV
			);
#endif

#if 0
		std::cout << std::endl;
		for (int i = 0; i < i_size; i++)
			std::cout << o_x[i] << ", ";
		std::cout << std::endl;
#endif
	}

};



#endif /* SRC_INCLUDE_LIBMATH_LAPACKBANDEDMATRIXSOLVER_HPP_ */
