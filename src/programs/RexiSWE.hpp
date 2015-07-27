/*
 * rexi_swe.hpp
 *
 *  Created on: 24 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_PROGRAMS_REXISWE_HPP_
#define SRC_PROGRAMS_REXISWE_HPP_

#if SWEET_USE_SPECTRAL_SPACE != 1
	#error	"Spectral space required for solvers"
#endif

#include <sweet/DataArray.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/SimulationParameters.hpp>
#include <complex>
#include <ostream>
#include <fstream>
#include <vector>

/**
 * This class implements the REXI (rational approximation of exponential integrator) solver for the SWE,
 * see High-order time-parallel approximation of evolution operators, T. Haut et. al.
 */
class RexiSWE
{
	typedef std::complex<double> complex;
	std::vector<complex> alpha_table;
	std::vector<complex> beta_table;

	double tau;
	int M;
	double f;


	void loadData(
			std::vector<complex> &i_data,
			bool i_load_real_values,
			const char *i_filename
	)
	{
		std::ifstream infile(i_filename);

		for (int i = 0; i < i_data.size(); i++)
		{
			if (!infile)
			{
				std::cerr << "Unexpected EOF in file " << i_filename << std::endl;
				exit(1);
			}

			double a;
			infile >> a;

			a /= 2.0*M_PIl;

			if (i_load_real_values)
				i_data[i].real(a);
			else
				i_data[i].imag(a);
		}
	}

	/**
	 * setup the REXI
	 */
public:
	void setup(
			double i_tau,	///< time step size
			int i_M,		///< number of sampling points
			double i_f		///< Coriolis force
	)
	{
		M = i_M;

		alpha_table.resize(M);
		beta_table.resize(M);

		std::ostringstream filename;

		filename << "data/rexi/alpha_" << M << "_re.txt";
		loadData(alpha_table, true, filename.str().c_str());
		filename .str("");
		filename .clear();

		filename << "data/rexi/alpha_" << M << "_im.txt";
		loadData(alpha_table, false, filename.str().c_str());
		filename .str("");
		filename .clear();

		filename << "data/rexi/beta_" << M << "_re.txt";
		loadData(beta_table, true, filename.str().c_str());
		filename .str("");
		filename .clear();

		filename << "data/rexi/beta_" << M << "_im.txt";
		loadData(beta_table, false, filename.str().c_str());
		filename .str("");
		filename .clear();

		tau = i_tau;

		f = i_f;
	}



	/**
	 * Solve the REXI of U(t) = exp(L*t)
	 *
	 * See also
	 * 	apply_rational_func_expL_wave_prop
	 * in Terrys code
	 */
	void run_timestep(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		Operators2D &op,
		const SimulationParameters &parameters
	)
	{
		if (parameters.sim_g != 1.0 || parameters.setup_h0 != 1.0)
		{
			std::cerr << "Only non-dimensional formulation supported yet, use g=1 and h0=1" << std::endl;
			exit(1);
		}

		DataArray<2> o_h(io_h.resolution);
		DataArray<2> o_u(io_h.resolution);
		DataArray<2> o_v(io_h.resolution);

		o_h.setAll(0);
		o_u.setAll(0);
		o_v.setAll(0);

#define TEST_IMAG_ZERO	1

#if TEST_IMAG_ZERO
		DataArray<2> o_h_im(o_h.resolution);
		DataArray<2> o_u_im(o_u.resolution);
		DataArray<2> o_v_im(o_v.resolution);

		o_h_im.setAll(0);
		o_u_im.setAll(0);
		o_v_im.setAll(0);
#endif

		for (int m = 0; m < M; m++)
		{
			// load alpha and scale by inverse of tau
			complex alpha = alpha_table[m]/tau;

			/*
			 * Note, that alpha and beta are complex valued!
			 *
			 * B = alpha^2 + f^2
			 *
			 *         ( alpha    -f  )
			 * A = 1/B (              )
			 *         (   f    alpha )
			 *
			 * D: (diff_x, diff_y) operator
			 *
			 * (D^2 - B) h = B/alpha ( h0 + D.(A.v0) )
			 */
			complex B = alpha*alpha + f*f;
			complex inv_B = 1.0/B;

			// complex derivative:
			// f = a+ib
			// df/dx = da/dx + i db/dx

			// from now on, we handle the complex and real valued parts separated.

			complex f_B = f/B;
			complex alpha_B = alpha/B;

			// assemble first element of A.v0
			DataArray<2> A_v0_l1_re = alpha_B.real()*io_v - f_B.real()*io_u;
			DataArray<2> A_v0_l1_im = alpha_B.imag()*io_v - f_B.imag()*io_u;

			// assemble second element of A.v0
			DataArray<2> A_v0_l2_re = f_B.real()*io_v + alpha_B.real()*io_u;
			DataArray<2> A_v0_l2_im = f_B.imag()*io_v + alpha_B.imag()*io_u;

			// compute
			// foo = h0 + D.(A.v0)
			DataArray<2> foo_re = parameters.setup_h0 + op.diff_c_x(A_v0_l1_re) + op.diff_c_y(A_v0_l2_re);
			DataArray<2> foo_im = op.diff_c_x(A_v0_l1_im) + op.diff_c_y(A_v0_l2_im);

			// bar = B/alpha
			complex bar = B/alpha;

			// compute rhs := foo*bar with complex numbers
			DataArray<2> rhs_re = foo_re*bar.real() - foo_im*bar.imag();
			DataArray<2> rhs_im = foo_re*bar.imag() + foo_im*bar.real();

			// compute
			// lhs/h := D^2 - B
			// split into real and imaginary parts
			DataArray<2> lhs_re = op.diff2_c_x + op.diff2_c_y - B.real();
			DataArray<2> lhs_im = op.diff2_c_x + op.diff2_c_y - B.imag();

			// now compute h = foobar / B
			DataArray<2> h_re(io_h.resolution);
			DataArray<2> h_im(io_h.resolution);

			{
				rhs_re.requestDataInCartesianSpace();
				rhs_im.requestDataInCartesianSpace();
				lhs_re.requestDataInCartesianSpace();
				lhs_im.requestDataInCartesianSpace();

				for (int j = 0; j < io_h.resolution[1]; j++)
				{
					for (int i = 0; i < io_h.resolution[0]; i++)
					{
						double ar = rhs_re.get(j, i);
						double ai = rhs_im.get(j, i);
						double br = lhs_re.get(j, i);
						double bi = lhs_im.get(j, i);

						double den = (br*br+bi*bi);

						if (den == 0)
						{
							// For Laplace operator, this is the integration constant C
							h_re.set(j, i, 0);
							h_im.set(j, i, 0);
						}
						else
						{
							h_re.set(j, i, (ar*br + ai*bi)/den);
							h_im.set(j, i, (ai*br - ar*bi)/den);
						}
					}
				}
			}


			// A.D h
			DataArray<2> Dx_h_re = op.diff_c_x(h_re);
			DataArray<2> Dy_h_re = op.diff_c_y(h_re);
			DataArray<2> Dx_h_im = op.diff_c_x(h_im);
			DataArray<2> Dy_h_im = op.diff_c_y(h_im);

			/*
			 *         ( alpha   -f   )
			 * A = 1/B |              |
			 *         (   f    alpha )
			 *
			 * v = −A.v0 + A.D h
			 */
			DataArray<2> u_re =
							- A_v0_l1_re
							+ alpha_B.real()*Dx_h_re - alpha_B.imag()*Dx_h_im
							- (f_B.real()*Dy_h_re - f_B.imag()*Dy_h_im)
						;
			DataArray<2> u_im =
							- A_v0_l1_im
							+ alpha_B.real()*Dx_h_im + alpha_B.real()*Dx_h_im
							- (f_B.real()*Dy_h_im + f_B.real()*Dy_h_im)
						;

			DataArray<2> v_re =
							- A_v0_l2_re
							+ f_B.real()*Dx_h_re - f_B.imag()*Dx_h_im
							- (alpha_B.real()*Dy_h_re - alpha_B.imag()*Dy_h_im)
						;
			DataArray<2> v_im =
							- A_v0_l2_im
							+ f_B.real()*Dx_h_im + f_B.real()*Dx_h_im
							- (alpha_B.real()*Dy_h_im + alpha_B.real()*Dy_h_im)
						;

			o_h += h_re*beta_table[m].real() - h_im*beta_table[m].imag();
			o_u += u_re*beta_table[m].real() - u_im*beta_table[m].imag();
			o_v += v_re*beta_table[m].real() - v_im*beta_table[m].imag();

#if TEST_IMAG_ZERO
			o_h_im += h_re*beta_table[m].imag() + h_im*beta_table[m].real();
			o_u_im += u_re*beta_table[m].imag() + u_im*beta_table[m].real();
			o_v_im += v_re*beta_table[m].imag() + u_im*beta_table[m].real();
#endif
		}

		double inv_tau = 1.0/tau;

		io_h = o_h * inv_tau;
		io_u = o_u * inv_tau;
		io_v = o_v * inv_tau;

		std::cout << "o_h : " << o_h.reduce_norm1_quad() << std::endl;
		std::cout << "o_u : " << o_u.reduce_norm1_quad() << std::endl;
		std::cout << "o_v : " << o_v.reduce_norm1_quad() << std::endl;

#if TEST_IMAG_ZERO
		o_h_im /= tau;
		o_u_im /= tau;
		o_v_im /= tau;

		std::cout << "o_h_im : " << o_h_im.reduce_norm1_quad() << std::endl;
		std::cout << "o_u_im : " << o_u_im.reduce_norm1_quad() << std::endl;
		std::cout << "o_v_im : " << o_v_im.reduce_norm1_quad() << std::endl;
#endif
	}

	~RexiSWE()
	{
	}
};

#endif /* SRC_PROGRAMS_REXISWE_HPP_ */
