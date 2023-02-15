/*
 * BandedMatrixSolver.hpp
 *
 *  Created on: 30 Jan 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 *              Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_LIBMATH_BANDEDMATRIXSOLVER_HPP_
#define SRC_INCLUDE_LIBMATH_BANDEDMATRIXSOLVER_HPP_

#include <complex>
#include <string.h>
#include <stdlib.h>
#include <limits>
#if SWEET_DEBUG
#include <chrono>
#endif

template <typename T>
class BandedMatrixSolverCommon
{
public:
	int max_N;
	int num_diagonals;
	int num_halo_size_diagonals;

	////int LDAB;

	///std::complex<double>* AB;
	int *IPIV;

	BandedMatrixSolverCommon()	:
		////AB(nullptr),
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
		this->max_N = i_max_N;
		this->num_diagonals = 2*i_num_off_diagonals+1;
		this->num_halo_size_diagonals = i_num_off_diagonals;

		assert(2*num_halo_size_diagonals+1 == num_diagonals);

		////LDAB = 2*num_halo_size_diagonals + num_halo_size_diagonals + 1;

		////AB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*LDAB*i_max_N);
		IPIV = (int*)malloc(sizeof(int)*i_max_N);
	}


	void shutdown()
	{
		if (IPIV != nullptr)
		{
			free(IPIV);
			IPIV = nullptr;
		}

		////if (AB != nullptr)
		////{
		////	free(AB);
		////	AB = nullptr;
		////}
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
class BandedMatrixSolver	:
		public BandedMatrixSolverCommon<T>
{
public:
	//////**
	///// * Solve a diagonal banded matrix
	///// *
	///// * A*X = B
	///// */
	/////void solve_diagBandedInverse_C(
	/////	const std::complex<double>* i_A,		///< Matrix for input and in-place transformations
	/////	const std::complex<double>* i_b,		///< RHS of equation and output of solution X
	/////	std::complex<double>* o_x,
	/////	int i_size
	/////);


	/////void solve_diagBandedInverse_Fortran(
	/////	const std::complex<double>* i_A,		///< Matrix for input and in-place transformations
	/////	const std::complex<double>* i_b,		///< RHS of equation and output of solution X
	/////	std::complex<double>* o_x,
	/////	int i_size
	/////);

	void solve_diagBandedInverse(
						int i_size,
						T* i_A,
						T* i_b,
						T* o_x,
						double i_small = 1e-15
			);



};


template <>
class BandedMatrixSolver<std::complex<double>>	:
		public BandedMatrixSolverCommon<std::complex<double>>
{
	typedef std::complex<double> T;

////////////////	/**
////////////////	 * Solve for input matrix
////////////////	 *
////////////////	 * i_A: cols: num_diagonals
////////////////	 *      rows: i_size
////////////////	 *
////////////////	 * The Fortran array will be transposed and with a size of (rows: LDAB, cols: i_size)
////////////////	 *
////////////////	 * i_b: RHS of equation
////////////////	 *
////////////////	 * o_x: Solution
////////////////	 */
////////////////public:
////////////////	void solve_diagBandedInverse_Carray(
////////////////		const std::complex<double>* i_A,
////////////////		const std::complex<double>* i_b,
////////////////		std::complex<double>* o_x,
////////////////		int i_size,
////////////////		int i_debug_block
////////////////	)	const
////////////////	{
////////////////		assert(max_N >= i_size);
////////////////
////////////////#if 0
////////////////		/**
////////////////		 * WARNING: LEAVE THIS BLOCK FOR DEBUGGING PURPOSE!!!
////////////////		 *
////////////////		 * 1) Decodes the compact stored matrix A into a full NxN matrix.
////////////////		 * 2) Convert to Fortran storage format
////////////////		 * 3) convert to LAPACK general band matrix format
////////////////		 */
////////////////		std::cout << "C ARRAY: A compactified" << std::endl;
////////////////		print_array_c(i_A, num_diagonals, i_size);
////////////////
////////////////		T *fortran_A = new T[i_size*i_size];
////////////////
////////////////		for (int i = 0; i < i_size*i_size; i++)
////////////////			fortran_A[i] = std::numeric_limits<double>::infinity();
////////////////
////////////////		// source c cols
////////////////		for (int i = 0; i < num_diagonals; i++)
////////////////		{
////////////////			// source c rows
////////////////			for (int j = 0; j < i_size; j++)
////////////////			{
////////////////				int di = j+i-num_halo_size_diagonals;
////////////////
////////////////				if (di < 0 || di >= i_size)
////////////////					continue;
////////////////
////////////////				int dj = j;
////////////////				fortran_A[di*i_size+dj] = i_A[j*num_diagonals+i];
////////////////			}
////////////////		}
////////////////
////////////////		std::cout << "FORTRAN ARRAY: A expanded" << std::endl;
////////////////		print_array_fortran(fortran_A, i_size, i_size);
////////////////
////////////////		for (int i = 0; i < i_size*LDAB; i++)
////////////////			AB[i] = std::numeric_limits<double>::infinity();
////////////////
////////////////		for (int i = 0; i < i_size; i++)
////////////////		{
////////////////			int fi = i+1;
////////////////
////////////////			// columns for output fortran array
////////////////			for (int j = 0; j < i_size; j++)
////////////////			{
////////////////				int fj = j+1;
////////////////
////////////////				if (std::max(1, fj-num_halo_size_diagonals) <= fi && fi <= std::min(i_size, fj+num_halo_size_diagonals))
////////////////					AB[(num_diagonals+fi-fj)-1 + (fj-1)*LDAB] = fortran_A[i+j*i_size];
////////////////			}
////////////////		}
////////////////
////////////////
////////////////		std::cout << "FORTRAN ARRAY: A compactified" << std::endl;
////////////////		print_array_fortran(AB, i_size, LDAB);
////////////////
////////////////
////////////////#if 1
////////////////		for (int i = 0; i < i_size*LDAB; i++)
////////////////			AB[i] = std::numeric_limits<double>::infinity();
////////////////
////////////////		// columns for output fortran array
////////////////		// rows for input c array
////////////////		for (int j = 0; j < i_size; j++)
////////////////		{
////////////////			// rows for output fortran array
////////////////			// columns for input c array
////////////////			for (int i = 0; i < num_diagonals; i++)
////////////////			{
////////////////				// compute square matrix indices
////////////////				int si = j+(num_halo_size_diagonals-i);
////////////////				int sj = j;
////////////////
//////////////////				std::cout << sj << " " << si << std::endl;
////////////////
////////////////				if (si < 0 || si >= i_size)
////////////////					continue;
////////////////
////////////////
////////////////				// AB is LDAB large!
////////////////				assert(LDAB*max_N > i*i_size+j);
////////////////				assert(LDAB*max_N > i+j*num_diagonals);
////////////////
////////////////				AB[(num_diagonals+si-sj-1) + sj*LDAB] = i_A[(j-i+num_halo_size_diagonals)*num_diagonals + i];
////////////////			}
////////////////		}
////////////////#endif
////////////////
////////////////		delete fortran_A;
////////////////
////////////////		std::cout << "FORTRAN ARRAY: A compactified (alternative, should match previous one)" << std::endl;
////////////////		print_array_fortran(AB, i_size, LDAB);
////////////////
////////////////#else
////////////////
////////////////#ifndef NDEBUG
////////////////		for (int i = 0; i < i_size*LDAB; i++)
////////////////			AB[i] = std::numeric_limits<double>::infinity();
////////////////#endif
////////////////
////////////////		// columns for output fortran array
////////////////		// rows for input c array
////////////////		for (int j = 0; j < i_size; j++)
////////////////		{
////////////////			// rows for output fortran array
////////////////			// columns for input c array
////////////////			for (int i = 0; i < num_diagonals; i++)
////////////////			{
////////////////				// compute square matrix indices
////////////////				int si = j+(num_halo_size_diagonals-i);
////////////////				int sj = j;
////////////////
////////////////				if (si < 0 || si >= i_size)
////////////////					continue;
////////////////
////////////////
////////////////				// AB is LDAB large!
////////////////				assert(LDAB*max_N > i*i_size+j);
////////////////				assert(LDAB*max_N > i+j*num_diagonals);
////////////////
////////////////				AB[(num_diagonals+si-sj-1) + sj*LDAB] = i_A[(j-i+num_halo_size_diagonals)*num_diagonals + i];
////////////////			}
////////////////		}
////////////////#endif
////////////////
////////////////		solve_diagBandedInverse_FortranArray(AB, i_b, o_x, i_size, i_debug_block);
////////////////	}



///////////////////////public:
///////////////////////	void solve_diagBandedInverse_FortranArray(
///////////////////////		const std::complex<double>* i_A,	///< A of max size
///////////////////////		const std::complex<double>* i_b,
///////////////////////		std::complex<double>* o_x,
///////////////////////		int i_size,
///////////////////////		int i_debug_block
///////////////////////	)	const
///////////////////////	{
///////////////////////		/*
///////////////////////		 * Make a copy of the array data since this is a destructive function
///////////////////////		 */
///////////////////////		if (AB != i_A)
///////////////////////			memcpy((void*)AB, (const void*)i_A, sizeof(std::complex<double>)*num_diagonals*LDAB);
///////////////////////
///////////////////////		memcpy((void*)o_x, (const void*)i_b, sizeof(std::complex<double>)*i_size);
///////////////////////
///////////////////////		solve_diagBandedInverse_FortranArray_inplace(AB, o_x, i_size, i_debug_block);
///////////////////////	}
///////////////////////
///////////////////////
///////////////////////
///////////////////////public:
///////////////////////	void solve_diagBandedInverse_FortranArray_inplace(
///////////////////////		std::complex<double>* io_A,		///< A of max size
///////////////////////		std::complex<double>* io_b_x,	///< rhs and solution x
///////////////////////		int i_size,
///////////////////////		int i_debug_block
///////////////////////	)	const
///////////////////////	{
///////////////////////		assert((num_diagonals & 1) == 1);
///////////////////////		assert(AB != nullptr);
///////////////////////
///////////////////////
///////////////////////////////////////#if 0
///////////////////////////////////////		std::cout << "************************************" << std::endl;
///////////////////////////////////////		std::cout << "i_size: " << i_size << std::endl;
///////////////////////////////////////		std::cout << "num_halo_size_diagonals: " << num_halo_size_diagonals << std::endl;
///////////////////////////////////////		std::cout << "LDAB: " << LDAB << std::endl;
///////////////////////////////////////#endif
///////////////////////////////////////
///////////////////////////////////////#if SWEET_LAPACK
///////////////////////////////////////		int info;
///////////////////////////////////////		zgbsv_(
///////////////////////////////////////				i_size,				// number of linear equations
///////////////////////////////////////				num_halo_size_diagonals,	// number of subdiagonals
///////////////////////////////////////				num_halo_size_diagonals,	// number of superdiagonals
///////////////////////////////////////				1,				// number of columns of matrix B
///////////////////////////////////////				io_A,				// array with matrix A to solve for
///////////////////////////////////////				LDAB,				// leading dimension of matrix A
///////////////////////////////////////				IPIV,				// integer array for pivoting
///////////////////////////////////////				io_b_x,				// output array
///////////////////////////////////////				i_size,				// leading dimension of array o_x
///////////////////////////////////////				info
///////////////////////////////////////			);
///////////////////////////////////////
///////////////////////////////////////		if (info != 0)
///////////////////////////////////////		{
///////////////////////////////////////			std::cerr << "Block ID: " << i_debug_block << std::endl;
///////////////////////////////////////			std::cerr << "zgbsv returned INFO != 0: " << info << std::endl;
///////////////////////////////////////			assert(false);
///////////////////////////////////////			exit(1);
///////////////////////////////////////		}
///////////////////////////////////////#else
///////////////////////////////////////		SWEETError("SWEET compiled without LAPACK!!!");
///////////////////////////////////////#endif
///////////////////////	}


public:

	int getMatrixIndex(
				int i_size,
				int i_i,
				int i_j
			) const
	{
		///return i_i * i_size + i_j;
		return (i_j - i_i + this->num_halo_size_diagonals) * this->num_diagonals + i_i;
	}

private:
	void swapRows(
			int i_size,
			T* io_A,
			T* io_b,
			int i_row1,
			int i_row2
		) const
	{
		for (int j = 0; j < i_size; j++)
			std::swap(io_A[getMatrixIndex(i_size, i_row1, j)], io_A[getMatrixIndex(i_size, i_row2, j)]);
		std::swap(io_b[i_row1], io_b[i_row2]);
	}

	void normalizeRow(
				int i_size,
				T* io_A,
				T* io_b,
				int* i_pivoting,
				int i_row,
				int i_ndiag_l = -1,
				int i_ndiag_u = -1
			) const
	{

		int col_min;
		int col_max;
		if (i_ndiag_l < 0)
			col_min = 0;
		else
			//col_min = std::max(0, i_pivoting[i_row] - i_ndiag_l);
			col_min = std::max(0, std::max(i_row, i_pivoting[i_row]) - i_ndiag_l);
		if (i_ndiag_u < 0)
			col_max = i_size;
		else
			//col_max = std::min(i_size, i_pivoting[i_row] + i_ndiag_u + 1);
		{
			// search highest i up to i_ndiag_l above
			int imax = i_pivoting[i_row];
			for (int i = 1; i < i_ndiag_l + 1; i++)
				if (i_row - i >= 0)
					if (i_pivoting[i_row - i] > imax)
						imax = i_pivoting[i_row - i];
			col_max = std::min(i_size, imax + i_ndiag_u + 1);
		}


		T fac = io_A[getMatrixIndex(i_size, i_row, i_row)];
		for (int j = col_min; j < col_max; j++)
			io_A[getMatrixIndex(i_size, i_row, j)] /= fac;
		io_b[i_row] /= fac;
	}

	void eliminateRowsBelow(
					int i_size,
					T* io_A,
					T* io_b,
					int* i_pivoting,
					int i_row,
					int i_ndiag_l = -1,
					int i_ndiag_u = -1
				) const
	{

		int row_max;
		if (i_ndiag_l < 0)
			row_max = i_size;
		else
			row_max = std::min(i_size, i_row + 1 + i_ndiag_l);
			/////row_max = i_size;

		for (int i = i_row + 1; i < row_max; i++)
		{

			int col_min;
			int col_max;
			if (i_ndiag_l < 0)
				col_min = i_row;
			else
				col_min = i_row;
				///col_min = std::max(i_row, std::min(i_pivoting[i_row], i_pivoting[i]) - i_ndiag_l);
			if (i_ndiag_u < 0)
				col_max = i_size;
			else
			{
				int imax = i_pivoting[i];
				for (int l = 1; l < i_ndiag_l + 1; l++)
					if (i - l >= 0)
						if (i_pivoting[i - l] > imax)
							imax = i_pivoting[i - l];

				col_max = std::min(i_size, imax + i_ndiag_u + 1);
				////col_max = i_size;
				////col_max = std::min(i_size, i + i_ndiag_u + 1);
				////std::cout << "AAAAAA " << i << " " << i_pivoting[i] << " " << col_min << " " << col_max << " " << imax << std::endl;
			}

			T fac = - io_A[getMatrixIndex(i_size, i, i_row)] / io_A[getMatrixIndex(i_size, i_row, i_row)];
			for (int j = col_min; j < col_max; j++)
				io_A[getMatrixIndex(i_size, i, j)] += fac * io_A[getMatrixIndex(i_size, i_row, j)];
			io_b[i] += fac * io_b[i_row];
		}

	}

	bool checkMatrixUpperTriangular(
						int i_size,
						T* i_A
					) const
	{
		bool is_upper_triangular = true;
		double small = 1e-15;

		for (int i = 0; i < i_size; i++)
		{
			for (int j = std::max(0, i - this->num_halo_size_diagonals); j < i; j++)
			{
				if (std::abs(i_A[getMatrixIndex(i_size, i , j)]) > small)
				{
					is_upper_triangular = false;
					break;
				}
			}
			if (!is_upper_triangular)
				break;
		}

		return is_upper_triangular;
	}

	void solveUpperTriangularSystem(
						int i_size,
						T* i_A,
						T* i_b,
						T* o_x
					) const
	{

		for (int i = i_size - 1; i >= 0; i--)
		{
			T rhs = i_b[i];
			for (int j = i + 1; j < i_size; j++)
				rhs -= o_x[j] * i_A[getMatrixIndex(i_size, i, j)];
			o_x[i] = rhs / i_A[getMatrixIndex(i_size, i, i)];
		}

	}

public:
	bool checkSolution(
					int i_size,
					T* i_A,
					T* i_b,
					T* i_x
				) const
	{
		bool solution_ok = true;
		double small = 1e-12;
		double max_diff = 0;

		for (int i = 0; i < i_size; i++)
		{
			T v = 0;
			for (int j = std::max(0, i - this->num_halo_size_diagonals); j < std::min(i_size, i + this->num_halo_size_diagonals + 1); j++)
			{
				std::cout << "AAAA " << i << " " << j << " " << i_x[j] << " " << i_A[getMatrixIndex(i_size, i, j)] << std::endl;
				v += i_x[j] * i_A[getMatrixIndex(i_size, i, j)];
			}
			max_diff = std::max(max_diff, std::abs(v - i_b[i]));
			std::cout << "v: " << v << std::endl;
			if (std::abs(v - i_b[i]) > small)
			{
				std::cout << std::setprecision(15) << "DIFF: " << std::abs(v - i_b[i]) << std::endl;
				solution_ok = false;
				break;
			}
		}

		std::cout << "MAX DIFF: " << max_diff << std::endl;

		if (!solution_ok)
		{
			this->printSolution(i_size, i_x);
		}

		return solution_ok;
	}

	void printAugmentedMatrix(
					int i_size,
					T* i_A,
					T* i_b
				) const
	{
		std::cout << "A | b" << std::endl;
		for (int i = 0; i < i_size; i++)
		{
			for (int j = 0; j < i_size; j++)
			{
				if (j >= std::max(0, i - this->num_halo_size_diagonals) && j <= std::min(i_size, i + this->num_halo_size_diagonals))
					std::cout << std::setprecision(3) << std::fixed << i_A[getMatrixIndex(i_size, i, j)] << "\t";
				else
					std::cout << std::setprecision(3) << std::fixed << (T)0. << "\t";
			}
			std::cout << "|\t" << i_b[i] << std::endl;
		}
	}

	void printSolution(
				int i_size,
				T* i_x
			) const
	{
		std::cout << "x = " << std::endl;
		for (int i = 0; i < i_size; i++)
			std::cout << i_x[i] << " ";
		std::cout << std::endl;
	}

	void copySystem(
				int i_size,
				const T* i_A,
				const T* i_b,
				T* &o_A2,
				T* &o_b2
			) const
	{
		for (int i = 0; i < i_size; i++)
		{
			o_b2[i] = i_b[i];
			for (int j = 0; j < i_size; j++)
			{
				int idx = getMatrixIndex(i_size, i, j);
				o_A2[idx] = i_A[idx];
			}
		}
	}

public:
	void solve_diagBandedInverse(
						int i_size,
						const T* i_A,
						const T* i_b,
						T* o_x,
#if SWEET_DEBUG
						double o_time_pivoting,
						double o_time_solving,
#endif
						double i_small = 1e-15
			) const
	{

		// copy A and b
		T* A = nullptr;
		T* b = nullptr;
		A = MemBlockAlloc::alloc<T>(i_size * this->num_diagonals * sizeof(T));
		b = MemBlockAlloc::alloc<T>(i_size * sizeof(T));
		this->copySystem(i_size, i_A, i_b, A, b);


#if SWEET_DEBUG
		////T* A2 = nullptr;
		////T* b2 = nullptr;
		////copySystem(i_size, i_A, i_b, A2, b2);

		auto time0 = std::chrono::high_resolution_clock::now();
#endif


		int ndiag_l = this->num_halo_size_diagonals;
		int ndiag_u = this->num_halo_size_diagonals;

		// loop over rows
		for (int i  = 0; i < i_size; i++)
		{
			// search for largest element in column
			int i_row_max = i;
			double val_max = std::abs(A[getMatrixIndex(i_size, i, i)]);
			bool pivoted_column = true;
			int lmax;
			if (ndiag_l < 0)
				lmax = i_size;
			else
				lmax = std::min(i_size, i + 1 + ndiag_l);
			for (int l = i + 1; l < lmax; l++)
			////for (int l = i + 1; l < i_size; l++)
			{
				int idx = getMatrixIndex(i_size, l, i);
				if (std::abs(A[idx]) > val_max)
				{
					val_max = std::abs(A[idx]);
					i_row_max = l;
				}
				if (std::abs(A[idx]) > i_small)
					pivoted_column = false;
			}

			// if all elements below are zero, nothing todo do
			if (pivoted_column)
				continue;

			// swap rows if needed
			if (i_row_max != i)
			{
				swapRows(i_size, A, b, i, i_row_max);
				std::swap(IPIV[i], IPIV[i_row_max]);
			}

			// normalize row
			normalizeRow(i_size, A, b, IPIV, i, ndiag_l, ndiag_u);

			// elimination in all rows below
			eliminateRowsBelow(i_size, A, b, IPIV, i, ndiag_l, ndiag_u);

		}
#if SWEET_DEBUG
		auto time1 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0);
		o_time_pivoting = duration.count();

		// check if matrix is upper triangular
		if (!this->checkMatrixUpperTriangular(i_size, A))
		{
			this->printAugmentedMatrix(i_size, A, b);
			SWEETError("Pivoted matrix is not upper triangular!");
		}
#endif

		///for (int i = 0; i < i_size; i++)
		///	std::cout << std::setprecision(15) << A[getMatrixIndex(i_size, i, i)] << std::endl;

		// solve upper triangular system
#if SWEET_DEBUG
		time0 = std::chrono::high_resolution_clock::now();
#endif

		this->solveUpperTriangularSystem(i_size, A, b, o_x);


		this->printAugmentedMatrix(i_size, A, b);

#if SWEET_DEBUG
		time1 = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0);
		o_time_solving = duration.count();
#endif

		delete [] A;
		delete [] b;

		MemBlockAlloc::free(A, i_size * this->num_diagonals * sizeof(T));
		MemBlockAlloc::free(b, i_size * sizeof(T));
	}

};



#endif /* SRC_INCLUDE_LIBMATH_BANDEDMATRIXSOLVER_HPP_ */
