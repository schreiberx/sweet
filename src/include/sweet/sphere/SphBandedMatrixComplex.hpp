/*
 * SPHSolverComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SPH_SPHSOLVER_COMPLEX_HPP_
#define SRC_INCLUDE_SPH_SPHSOLVER_COMPLEX_HPP_

#include "../libmath/LapackBandedMatrixSolver.hpp"
#include "../sphere/SphereSPHIdentities.hpp"
#include "BandedMatrixComplex.hpp"



/**
 * phi(lambda,mu) denotes the solution
 */
template <typename T>	// T: complex valued single or double precision method
class SphBandedMatrixComplex	:
		SphereSPHIdentities
{
public:
	/**
	 * Matrix on left-hand side
	 */
	BandedMatrixComplex<T> lhs;

	/**
	 * SPH configuration
	 */
	SphereDataConfig *sphConfig;

	/**
	 * Solver for banded matrix
	 */
	SphBandedMatrix< std::complex<double> > bandedMatrixSolver;

	/**
	 * Buffer to compactify the N-varying for a particular M
	 */
	std::complex<double> *buffer_in, *buffer_out;


	/**
	 * Setup the SPH solver
	 */
public:
	void setup(
			SphereDataConfig *i_sphConfig,		///< Handler to sphConfig
			int i_halosize_offdiagonal	///< Size of the halo around. A value of 2 allocates data for 5 diagonals.
	)
	{
		sphConfig = i_sphConfig;

		lhs.setup(sphConfig, i_halosize_offdiagonal);

		bandedMatrixSolver.setup(i_sphConfig->spec_n_max+1, i_halosize_offdiagonal);

		buffer_in = MemBlockAlloc::alloc< std::complex<double> >(sphConfig->spec_n_max+1);
		buffer_out = MemBlockAlloc::alloc< std::complex<double> >(sphConfig->spec_n_max+1);
	}


	SphBandedMatrixComplex()	:
		sphConfig(nullptr),
		buffer_in(nullptr),
		buffer_out(nullptr)
	{
	}


	~SphBandedMatrixComplex()
	{
		if (buffer_in != nullptr)
		{
			MemBlockAlloc::free(buffer_in, sphConfig->spec_num_elems);
			MemBlockAlloc::free(buffer_out, sphConfig->spec_num_elems);
		}
	}


	/**
	 * Solver for
	 * 	a*phi(lambda,mu)
	 */
	void solver_component_scalar_phi(
			const std::complex<double> &i_value
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, 0, i_value);
			}
		}
	}



	/**
	 * Solver for
	 * 	mu*phi(lambda,mu)
	 */
	void solver_component_mu_phi(
			const std::complex<double> &i_scalar = 1.0
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -1, R(n-1,m)*i_scalar);
				lhs.rowElement_add(row, n, m, +1, S(n+1,m)*i_scalar);
			}
		}
	}



	/**
	 * Solver for
	 * 	(1-mu*mu)*d/dmu phi(lambda,mu)
	 */
	void solver_component_one_minus_mu_mu_diff_mu_phi()
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -1, (-(double)n+1.0)*R(n-1,m));
				lhs.rowElement_add(row, n, m, +1, ((double)n+2.0)*S(n+1,m));
			}
		}
	}



	/**
	 * Solver for
	 * 	scalar*phi(lambda,mu)
	 */
	void solver_component_rexi_z1(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		solver_component_scalar_phi(i_scalar);
	}



	/**
	 * Solver for
	 * 	mu^2*phi(lambda,mu)
	 */
	void solver_component_rexi_z2(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
#if 0
				lhs.rowElement_add(row, n, m, -2, i_scalar*R(n-1,m)*R(n-2,m));
				lhs.rowElement_add(row, n, m,  0, i_scalar*R(n-1,m)*S(n,m) + S(n+1,m)*R(n,m));
				lhs.rowElement_add(row, n, m, +2, i_scalar*S(n+1,m)*S(n+2,m));
#else
				lhs.rowElement_add(row, n, m, -2, i_scalar*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, i_scalar*B(n,m));
				lhs.rowElement_add(row, n, m, +2, i_scalar*C(n+2,m));
#endif
			}
		}
	}



	/**
	 * Solver for
	 * 	mu^4*phi(lambda,mu)
	 */
	void solver_component_rexi_z3(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m, -4, i_scalar*(A(n-2,m)*A(n-4.0,m)));
				lhs.rowElement_add(row, n, m, -2, i_scalar*(A(n-2,m)*B(n-2,m)+B(n,m)*A(n-2,m)));
				lhs.rowElement_add(row, n, m,  0, i_scalar*(A(n-2,m)*C(n,m)+B(n,m)*B(n,m)+C(n+2,m)*A(n,m)));
				lhs.rowElement_add(row, n, m, +2, i_scalar*(B(n,m)*C(n+2,m)+C(n+2,m)*B(n+2,m)));
				lhs.rowElement_add(row, n, m, +4, i_scalar*(C(n+2,m)*C(n+4,m)));
			}
		}
	}


	/**
	 * Solver for
	 * Z4 := grad_j(mu) * grad_i(phi)
	 */
	void solver_component_rexi_z4(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  0, 1.0/(i_r*i_r)*i_scalar*T(0, m));
			}
		}
	}



	/**
	 * Solver for
	 * Z4robert := grad_j(mu) * grad_i(phi)
	 */
	void solver_component_rexi_z4robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			std::complex<double> fac = 1.0/(i_r*i_r)*i_scalar*T(0, m);
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -2, -fac*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, -fac*B(n,m) + fac);
				lhs.rowElement_add(row, n, m, +2, -fac*C(n+2,m));
			}
		}
	}



	/**
	 * Solver for
	 * Z5 := grad_j(mu) * mu^2 * grad_i(phi)
	 */
	void solver_component_rexi_z5(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				std::complex<double> fac = 1.0/(i_r*i_r)*i_scalar*T(0, m);

				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  -2, fac*A(n-2,m)	);
				lhs.rowElement_add(row, n, m,   0, fac*B(n+0,m)	);
				lhs.rowElement_add(row, n, m,  +2, fac*C(n+2,m)	);
			}
		}
	}


	/**
	 * Solver for
	 * Z5robert := grad_j(mu) * mu^2 * grad_i(phi)
	 */
	void solver_component_rexi_z5robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			std::complex<double> fac = 1.0/(i_r*i_r)*i_scalar*T(0, m);
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  -2,    -fac*(A(n-2, m)*A(n-4, m))    );
				lhs.rowElement_add(row, n, m,  -2,    -fac*(A(n-2, m)*B(n-2, m) + B(n, m)*A(n-2, m))                 + fac*A(n-2, m)        );
				lhs.rowElement_add(row, n, m,  -2,    -fac*(A(n-2, m)*C(n, m) + B(n, m)*B(n, m) + C(n+2, m)*A(n, m)) + fac*B(n, m)          );
				lhs.rowElement_add(row, n, m,  -2,    -fac*(B(n, m)*C(n+2, m) + C(n+2, m)*B(n+2, m))                 + fac*C(n+2, m)        );
				lhs.rowElement_add(row, n, m,  -2,    -fac*(C(n+2, m)*C(n+4, m))    );
			}
		}
	}



	/**
	 * Solver for
	 * Z6 := grad_j(mu) * mu * grad_j(phi)
	 */
	void solver_component_rexi_z6(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  -2, -1.0/(i_r*i_r)*i_scalar*(D(n-1,m)*R(n-2,m) + (E(n,m)-3.0)*A(n-2,m)));
				lhs.rowElement_add(row, n, m,   0, -1.0/(i_r*i_r)*i_scalar*(1.0+D(n-1,m)*S(n,m)+(E(n,m)-3.0)*B(n,m)));
				lhs.rowElement_add(row, n, m,  +2, -1.0/(i_r*i_r)*i_scalar*(E(n,m)-3.0)*C(n+2,m));
			}
		}
	}

	/**
	 * Solver for
	 * Z6robert := grad_j(mu) * mu * grad_j(phi)
	 */
	void solver_component_rexi_z6robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			std::complex<double> s = 1.0/(i_r*i_r)*i_scalar*(double)m;
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m,  -2, s*(G(n-1,m)*R(n-2,m))	);
				lhs.rowElement_add(row, n, m,   0, s*(G(n-1,m)*S(n,m) + H(n+1,m)*R(n,m))	);
				lhs.rowElement_add(row, n, m,  +2, s*(H(n+1,m)*S(n+2,m))	);

				lhs.rowElement_add(row, n, m,  -4, s*(-A(n-4,m)*(G(n-1,m)*R(n-2,m)))	);
				lhs.rowElement_add(row, n, m,  -2, s*(-B(n-2,m)*(G(n-1,m)*R(n-2,m)))	);
				lhs.rowElement_add(row, n, m,   0, s*(-C(n+0,m)*(G(n-1,m)*R(n-2,m)))	);

				lhs.rowElement_add(row, n, m,  -2, s*(-A(n-2,m)*(G(n-1,m)*S(n,m) + H(n+1,m)*R(n,m)))	);
				lhs.rowElement_add(row, n, m,   0, s*(-B(n+0,m)*(G(n-1,m)*S(n,m) + H(n+1,m)*R(n,m)))	);
				lhs.rowElement_add(row, n, m,  +2, s*(-C(n+2,m)*(G(n-1,m)*S(n,m) + H(n+1,m)*R(n,m)))	);

				lhs.rowElement_add(row, n, m,   0, s*(-A(n+0,m)*(H(n+1,m)*S(n+2,m)))	);
				lhs.rowElement_add(row, n, m,  +2, s*(-B(n+2,m)*(H(n+1,m)*S(n+2,m)))	);
				lhs.rowElement_add(row, n, m,  +4, s*(-C(n+4,m)*(H(n+1,m)*S(n+2,m)))	);
			}
		}
	}


	/**
	 * Solver for
	 * Z7 := laplace(phi)
	 */
	void solver_component_rexi_z7(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, 0, -1.0/(i_r*i_r)*i_scalar*(double)n*((double)n+1.0));
			}
		}
	}


	/**
	 * Solver for
	 * Z8 := mu*mu*laplace(phi)
	 */
	void solver_component_rexi_z8(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#pragma omp parallel for
		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m, 0, 1.0/(i_r*i_r)*i_scalar*2.0);

				double fac = -1.0/(i_r*i_r)*((double)m*m-(double)n*(n+1.0));
				lhs.rowElement_add(row, n, m, -2, i_scalar*fac*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, i_scalar*fac*B(n,m));
				lhs.rowElement_add(row, n, m, +2, i_scalar*fac*C(n+2,m));
			}
		}
	}



	/**
	 * Apply the solver matrix.
	 * This function is intended to be used for debugging.
	 * WARNING: This only multiplies the i_x values with the matrix.
	 * Use solve(...) to solve for the matrix
	 */
	SphereDataComplex apply(
			const SphereDataComplex &i_x	///< solution to be searched
	)
	{
		SphereDataComplex out(sphConfig);

#pragma omp parallel for
		for (int n = 0; n <= sphConfig->spec_n_max; n++)
		{
			int idx = sphConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out.data_spec[idx] = 0;

				std::complex<double> *row = lhs.getMatrixRow(n, m);
				for (int i = 0; i < lhs.num_diagonals; i++)
				{
					int delta = i-lhs.halosize_off_diagonal;
					out.data_spec[idx] += lhs.rowElement_getRef(row, n, m, delta)*i_x.spec_get(n+delta, m);
				}

				idx++;
			}
		}

		out.data_spat_valid = false;
		out.data_spec_valid = true;

		return out;
	}


	SphereDataComplex solve(
			const SphereDataComplex &i_rhs
	)
	{
		i_rhs.request_data_spectral();

		SphereDataComplex out(sphConfig);

		for (int m = -sphConfig->spec_m_max; m <= sphConfig->spec_m_max; m++)
		{
			int idx = sphConfig->getArrayIndexByModes_Complex_NCompact(std::abs(m),m);

			/*
			 * Compactify RHS coefficients
			 */
			{
				int buffer_idx = 0;
				for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
				{
					buffer_in[buffer_idx] = i_rhs.data_spec[sphConfig->getArrayIndexByModes_Complex(n,m)];
					buffer_idx++;
				}
			}

			bandedMatrixSolver.solve_diagBandedInverse_Carray(
							&lhs.data[idx*lhs.num_diagonals],
							buffer_in,
							buffer_out,
							sphConfig->spec_n_max+1-std::abs(m)	// size of block (same as for SPHSolver)
					);


			/*
			 * Compactify RHS coefficients
			 */
			{
				int buffer_idx = 0;
				for (int n = std::abs(m); n <= sphConfig->spec_n_max; n++)
				{
					out.data_spec[sphConfig->getArrayIndexByModes_Complex(n,m)] = buffer_out[buffer_idx];
					buffer_idx++;
				}
			}
		}

		out.data_spec_valid = true;
		out.data_spat_valid = false;

		return out;
	}
};


#endif /* SRC_INCLUDE_SPH_SPHSOLVER_HPP_ */
