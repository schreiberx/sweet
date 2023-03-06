/*
 * SPHSolverComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SPH_BANDED_MATRIX_COMPLEX_HPP_
#define SRC_INCLUDE_SPH_BANDED_MATRIX_COMPLEX_HPP_

#include <sweet/libmath/BandedMatrixPhysicalComplex.hpp>
#include <sweet/libmath/LapackBandedMatrixSolver.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/core/sphere/SphereHelpers_SPHIdentities.hpp>



/**
 * phi(lambda,mu) denotes the solution
 */
template < typename T = std::complex<double> >	// T: complex valued single or double precision method
class SphBandedMatrixPhysicalComplex	:
		sweet::SphereHelpers_SPHIdentities
{
public:
	/**
	 * Matrix on left-hand side
	 */
	BandedMatrixPhysicalComplex<T> lhs;

	/**
	 * SPH configuration
	 */
	const sweet::SphereDataConfig *sphereDataConfig;

	/**
	 * Solver for banded matrix
	 */
	LapackBandedMatrixSolver< std::complex<double> > bandedMatrixSolver;

	/**
	 * Size of buffers
	 */
	std::size_t buffer_size;

	/**
	 * Buffer to compactify the N-varying for a particular M
	 */
	std::complex<double> *buffer_in, *buffer_out;

	/**
	 * Setup the SPH solver
	 */
public:
	void setup(
			const sweet::SphereDataConfig *i_sphereDataConfig,		///< Handler to sphereDataConfig
			int i_halosize_offdiagonal	///< Size of the halo around. A value of 2 allocates data for 5 diagonals.
	)
	{
		sphereDataConfig = i_sphereDataConfig;

		lhs.setup(sphereDataConfig, i_halosize_offdiagonal);

		bandedMatrixSolver.setup(i_sphereDataConfig->spectral_modes_n_max+1, i_halosize_offdiagonal);

		buffer_size = (sphereDataConfig->spectral_modes_n_max+1)*sizeof(std::complex<double>);

		buffer_in = MemBlockAlloc::alloc< std::complex<double> >(buffer_size);
		buffer_out = MemBlockAlloc::alloc< std::complex<double> >(buffer_size);
	}


	SphBandedMatrixPhysicalComplex()	:
		sphereDataConfig(nullptr),
		buffer_size(0),
		buffer_in(nullptr),
		buffer_out(nullptr)
	{
	}


	~SphBandedMatrixPhysicalComplex()
	{
		if (buffer_in != nullptr)
		{
			MemBlockAlloc::free(buffer_in, buffer_size);
			MemBlockAlloc::free(buffer_out, buffer_size);
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
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
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
#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
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
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
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
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m, -2, i_scalar*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, i_scalar*B(n,m));
				lhs.rowElement_add(row, n, m, +2, i_scalar*C(n+2,m));
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
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m, -4, i_scalar*(A(n-2,m)*A(n-4,m)));
				lhs.rowElement_add(row, n, m, -2, i_scalar*(A(n-2,m)*B(n-2,m)+B(n  ,m)*A(n-2,m)));
				lhs.rowElement_add(row, n, m,  0, i_scalar*(A(n-2,m)*C(n,  m)+B(n  ,m)*B(n  ,m)+C(n+2,m)*A(n,m)));
				lhs.rowElement_add(row, n, m, +2, i_scalar*(B(n,  m)*C(n+2,m)+C(n+2,m)*B(n+2,m)));
				lhs.rowElement_add(row, n, m, +4, i_scalar*(C(n+2,m)*C(n+4,m)));
			}
		}
	}


#if 0
	/**
	 * Solver for
	 * 	mu^4*phi(lambda,mu)
	 */
	void solver_component_rexi_z3c(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
//			std::cout << "**********************************" << std::endl;
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

//				std::cout << std::endl;
//				std::cout << Ac(n-2,m) << "\t" << Ac(n-4,m) << std::endl;
				lhs.rowElement_add(row, n, m, -4, i_scalar*(Ac(n-2,m)*Ac(n-4,m)));

//				std::cout << Ac(n-2,m)*Bc(n-2,m) << "\t" << Bc(n  ,m)*Ac(n-2,m) << std::endl;
				lhs.rowElement_add(row, n, m, -2, i_scalar*(Ac(n-2,m)*Bc(n-2,m) + Bc(n  ,m)*Ac(n-2,m)));

//				std::cout << Ac(n-2,m)*Cc(n,  m) << "\t" << Bc(n  ,m)*Bc(n  ,m) << "\t" << Cc(n+2,m)*Ac(n,m) << std::endl;
				lhs.rowElement_add(row, n, m,  0, i_scalar*(Ac(n-2,m)*Cc(n,  m) + Bc(n  ,m)*Bc(n  ,m) + Cc(n+2,m)*Ac(n,m)));

//				std::cout << Bc(n,  m)*Cc(n+2,m) << "\t" << Cc(n+2,m)*Bc(n+2,m) << std::endl;
				lhs.rowElement_add(row, n, m, +2, i_scalar*(Bc(n,  m)*Cc(n+2,m) + Cc(n+2,m)*Bc(n+2,m)));

//				std::cout << Cc(n+2,m)*Cc(n+4,m) << std::endl;
				lhs.rowElement_add(row, n, m, +4, i_scalar*(Cc(n+2,m)*Cc(n+4,m)));
			}
		}
	}
#endif


	/**
	 * Solver for
	 * Z4 := grad_j(mu) * grad_i(phi)
	 */
	void solver_component_rexi_z4(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  0, 1.0/(i_r*i_r)*i_scalar*T(0, m));
			}
		}
	}



	/**
	 * Solver for
	 * Z4robert := grad_j(mu) * grad_i(phi)
	 *           = im/(r*r) * F_n^m(Phi)
	 */
	void solver_component_rexi_z4robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::complex<double> fac = 1.0/(i_r*i_r)*i_scalar*T(0, m);
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  0, fac);
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
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			T fac = std::complex<double>(1.0/(i_r*i_r))*i_scalar*T(0, m);

			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  -2, fac*A(n-2,m)	);
				lhs.rowElement_add(row, n, m,   0, fac*B(n+0,m)	);
				lhs.rowElement_add(row, n, m,  +2, fac*C(n+2,m)	);
			}
		}
	}


	/**
	 * Solver for
	 *
	 * Z5robert := grad_j(mu) * mu^2 * grad_i(phi)
	 *           = i*m/(r*r) * F_n^m(mu^2 \phi)
	 */
	void solver_component_rexi_z5robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::complex<double> fac = (i_scalar/(i_r*i_r))*std::complex<double>(0, m);

			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -2, fac*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, fac*B(n,m));
				lhs.rowElement_add(row, n, m, +2, fac*C(n+2,m));
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
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
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
	 *           =
	 */
	void solver_component_rexi_z6robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		/*
		 * First part
		 */
		// phi
		solver_component_rexi_z1(-i_scalar/(i_r*i_r), i_r);

		// mu^2*phi
		solver_component_rexi_z2(3.0*i_scalar/(i_r*i_r), i_r);


		/*
		 * Second part
		 */
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			std::complex<double> fac = -1.0/(i_r*i_r)*i_scalar;
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m,  -2, fac*(G(n-1,m)*R(n-2,m))	);
				lhs.rowElement_add(row, n, m,   0, fac*(G(n-1,m)*S(n,m) + H(n+1,m)*R(n,m))	);
				lhs.rowElement_add(row, n, m,  +2, fac*(H(n+1,m)*S(n+2,m))	);
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
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
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
		std::complex<double> fac = (1.0/(i_r*i_r))*i_scalar;

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// eq (D) in REXI SH document section 8.1.11 ver. 15
				std::complex<double> a = fac*(-(double)n*((double)n+1.0));
				lhs.rowElement_add(row, n, m, -2, a*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, a*B(n,m));
				lhs.rowElement_add(row, n, m, +2, a*C(n+2,m));

				// eq (Ba+Aa) in REXI SH document section 8.1.11 ver. 15
				lhs.rowElement_add(row, n, m, -2, 4.0*fac*( R(n-2,m)*G(n-1,m)) );
				lhs.rowElement_add(row, n, m,  0, 4.0*fac*( S(n,m)*G(n-1,m) + R(n,m)*H(n+1,m)) );
				lhs.rowElement_add(row, n, m, +2, 4.0*fac*( S(n+2,m)*H(n+1,m)) );

				// eq (Ab) in REXI SH document section 8.1.11 ver. 15
				lhs.rowElement_add(row, n, m, 0, 2.0*fac);

				lhs.rowElement_add(row, n, m, -2, -6.0*fac*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, -6.0*fac*B(n,m));
				lhs.rowElement_add(row, n, m, +2, -6.0*fac*C(n+2,m));
			}
		}
	}




	/**
	 * Solver for implicit time integration
	 */
	void solver_component_implicit_FJinvF(
			const std::complex<double>& i_dt_two_omega
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// out of boundary check for P(n-2, m)
				if (n-2 >= std::abs(m))
				{
					Tcomplex a =
							i_dt_two_omega * i_dt_two_omega
							* implicit_f_minus(n, m)
							/ implicit_J_scalar(n-1, m, i_dt_two_omega)
							* implicit_f_minus(n-1, m);

					lhs.rowElement_add_NEW(row, n, m, -2, a);
				}

				{
					// out of boundary check for P(n, m) (at least something like that)
					if (n > std::abs(m))
					{
						Tcomplex a =
								i_dt_two_omega * i_dt_two_omega
								* implicit_f_minus(n, m)
								/ implicit_J_scalar(n-1, m, i_dt_two_omega)
								* implicit_f_plus(n-1, m);

						lhs.rowElement_add_NEW(row, n, m, 0, a);
					}

					// out of boundary check for P(n, m) (at least something like that)
					if (n < sphereDataConfig->spectral_modes_n_max)
					{
						Tcomplex a =
								i_dt_two_omega * i_dt_two_omega
								* implicit_f_plus(n, m)
								/ implicit_J_scalar(n+1, m, i_dt_two_omega)
								* implicit_f_minus(n+1, m);

						lhs.rowElement_add_NEW(row, n, m, 0, a);
					}

				}

				// out of boundary check for P(n+2, m)
				if (n+2 <= sphereDataConfig->spectral_modes_n_max)
				{
					Tcomplex a =
							i_dt_two_omega * i_dt_two_omega
							* implicit_f_plus(n, m)
							/ implicit_J_scalar(n+1, m, i_dt_two_omega)
							* implicit_f_plus(n+1, m);

					lhs.rowElement_add_NEW(row, n, m, 2, a);
				}
			}
		}
	}



	/**
	 * Solver for implicit time integration
	 */
	void solver_component_implicit_F(
			const std::complex<double>& i_dt_two_imega
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// out of boundary check for P(n-1, m)
				if (n-1 >= std::abs(m))
				{
					Tcomplex a = implicit_f_minus(n, m) * i_dt_two_imega;
					lhs.rowElement_add_NEW(row, n, m, -1, a);
				}

				// out of boundary check for P(n+1, m)
				if (n+1 <= sphereDataConfig->spectral_modes_n_max)
				{
					Tcomplex a = implicit_f_plus(n, m) * i_dt_two_imega;

					lhs.rowElement_add_NEW(row, n, m, 1, a);
				}
			}
		}
	}



	/**
	 * Solver for implicit time integration
	 */
	void solver_component_implicit_FJinv(
			const std::complex<double>& i_dt_two_imega
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// out of boundary check for P(n-1, m)
				if (n-1 >= std::abs(m))
				{
					Tcomplex a =	i_dt_two_imega
									* implicit_f_minus(n, m)
									/ implicit_J_scalar(n-1, m, i_dt_two_imega);

					lhs.rowElement_add_NEW(row, n, m, -1, a);
				}

				// out of boundary check for P(n+1, m)
				if (n+1 <= sphereDataConfig->spectral_modes_n_max)
				{
					// sort out div/0 from implicit_J_scalar
					//if (n != 0)
					{
						Tcomplex a =	i_dt_two_imega
										* implicit_f_plus(n, m)
										/ implicit_J_scalar(n+1, m, i_dt_two_imega);

						lhs.rowElement_add_NEW(row, n, m, 1, a);
					}
				}
			}
		}
	}



	/**
	 * Solver for implicit time integration
	 */
	void solver_component_implicit_J(
			const std::complex<double>& i_dt_two_imega
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				Tcomplex a = implicit_J_scalar(n, m, i_dt_two_imega);
				lhs.rowElement_add_NEW(row, n, m, 0, a);
			}
		}
	}



	/**
	 * Solver for implicit time integration
	 */
	void solver_component_implicit_I(
			const std::complex<double>& i_scalar
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add_NEW(row, n, m, 0, i_scalar);
			}
		}
	}



	/**
	 * Solver for implicit time integration
	 */
	void solver_component_implicit_L(
			const std::complex<double>& i_scalar,
			const std::complex<double>& i_dt,		// note, that the dt is also in L. That's not a bug
			const double i_radius
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				Tcomplex a = i_scalar*i_dt*(n*(n+1.0))/(i_radius*i_radius);
				lhs.rowElement_add_NEW(row, n, m, 0, a);
			}
		}
	}





	/**
	 * Apply the solver matrix.
	 * This function is intended to be used for debugging.
	 * WARNING: This only multiplies the i_x values with the matrix.
	 * Use solve(...) to solve for the matrix
	 */
	sweet::SphereData_SpectralComplex apply(
			const sweet::SphereData_SpectralComplex &i_x	///< solution to be searched
	)
	{
		sweet::SphereData_SpectralComplex out(sphereDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out.spectral_space_data[idx] = 0;

				std::complex<double> *row = lhs.getMatrixRow(n, m);
				for (int i = 0; i < lhs.num_diagonals; i++)
				{
					int delta = i-lhs.halosize_off_diagonal;
					out.spectral_space_data[idx] += lhs.rowElement_getRef(row, n, m, delta)*i_x.spectral_get_(n+delta, m);
				}

				idx++;
			}
		}

		return out;
	}




	sweet::SphereData_SpectralComplex solve(
			const sweet::SphereData_PhysicalComplex &i_rhs
	)
	{
		sweet::SphereData_SpectralComplex tmp = i_rhs;
		return solve(tmp);
	}

	sweet::SphereData_SpectralComplex solve(
			const sweet::SphereData_SpectralComplex &i_rhs
	)
	{
		sweet::SphereData_SpectralComplex out(sphereDataConfig);

		i_rhs.check_sphereDataConfig_identical_res(sphereDataConfig);

		for (int m = -sphereDataConfig->spectral_modes_m_max; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			int idx = sphereDataConfig->getArrayIndexByModes_Complex_NCompact(std::abs(m),m);

			/*
			 * Compactify RHS coefficients
			 */
			{
				int buffer_idx = 0;
				for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
				{
					assert(buffer_idx < sphereDataConfig->spectral_modes_n_max+1);
					buffer_in[buffer_idx] = i_rhs.spectral_space_data[sphereDataConfig->getArrayIndexByModes_Complex(n,m)];
					buffer_idx++;
				}
			}

			bandedMatrixSolver.solve_diagBandedInverse_Carray(
							&lhs.data[idx*lhs.num_diagonals],
							buffer_in,
							buffer_out,
							sphereDataConfig->spectral_modes_n_max+1-std::abs(m),	// size of block (same as for SPHSolver)
							idx
					);


			/*
			 * Compactify RHS coefficients
			 */
			{
				int buffer_idx = 0;
				for (int n = std::abs(m); n <= sphereDataConfig->spectral_modes_n_max; n++)
				{
					out.spectral_space_data[sphereDataConfig->getArrayIndexByModes_Complex(n,m)] = buffer_out[buffer_idx];
					buffer_idx++;
				}
			}
		}

		return out;
	}
};


#endif
