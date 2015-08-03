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

	std::vector<complex> filter_alpha_table;
	std::vector<complex> filter_beta_table;

	double tau;
	int M;
	int M_filter;
	double f;


	void loadData(
			std::vector<complex> &i_data,
			bool i_load_real_values,
			const char *i_filename
	)
	{
		std::ifstream infile(i_filename);

		for (std::size_t i = 0; i < i_data.size(); i++)
		{
			if (!infile)
			{
				std::cerr << "Unexpected EOF in file " << i_filename << std::endl;
				exit(1);
			}

			double a;
			infile >> a;

//			a /= 2.0*M_PIl;

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
			int i_filter_M,	///< number of points for filtering
			double i_f		///< Coriolis force
	)
	{
		M = i_M;

		alpha_table.resize(M);
		beta_table.resize(M);

		std::ostringstream filename;

		filename << "data/rexi/rexi_" << M << "_poles_re.txt";
		loadData(alpha_table, true, filename.str().c_str());
		filename .str("");
		filename .clear();

		filename << "data/rexi/rexi_" << M << "_poles_im.txt";
		loadData(alpha_table, false, filename.str().c_str());
		filename .str("");
		filename .clear();

		filename << "data/rexi/rexi_" << M << "_coef_re.txt";
		loadData(beta_table, true, filename.str().c_str());
		filename .str("");
		filename .clear();

		filename << "data/rexi/rexi_" << M << "_coef_im.txt";
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
	 *
	 * See also doc/rexi/understanding_rexi.pdf
	 */
	void run_timestep(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		Operators2D &op,
		const SimulationParameters &i_parameters
	)
	{
		if (i_parameters.sim_g != 1.0 || i_parameters.setup_h0 != 1.0)
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
			/*
			 * (D^2 + K) h = k/a h(0) - H f/a D x v(0) + H D.v(0)
			 */
			// load alpha and scale by inverse of tau
			complex alpha = alpha_table[m]/tau;
			complex kappa = alpha*alpha + f*f;

			// TERM 1: k/a h(0)
			complex k_a = kappa/alpha;

			DataArray<2> rhs_re = k_a.real() * io_h;
			DataArray<2> rhs_im = k_a.imag() * io_h;

			// TERM 2: H f/a D x v(0)
			complex H_f_a = i_parameters.setup_h0*i_parameters.sim_f/alpha;
			DataArray<2> Dxv = op.diff_c_x(io_v) - op.diff_c_y(io_u);
			rhs_re -= H_f_a.real()*Dxv;
			rhs_im -= H_f_a.imag()*Dxv;

			// TERM 3: H D.v(0)
			DataArray<2> Dv = op.diff_c_x(io_u) + op.diff_c_y(io_v);
			rhs_re += i_parameters.setup_h0*Dv;

			// Now we have the RHS assembled.
			// It's stored in Cartesian space and is split into imaginary and real components.
			// In order to be able to compute the inverse of (D^2-k), we first have to compute
			// the FFT on these complex numbers.

			Complex2DArrayFFT rhs_complex(rhs_re.resolution);
			rhs_complex.loadRealAndImagFromDataArrays(rhs_re, rhs_im);

			// assemble LHS (D^2-k)
			// lhs/h := (D^2 - B)*h
			Complex2DArrayFFT lhs_complex(rhs_re.resolution);
			lhs_complex.loadRealFromDataArray(op.diff2_c_x + op.diff2_c_y - kappa.real());
			lhs_complex.setAllIm(-kappa.imag());

			// compute h:=RHS/LHS
			Complex2DArrayFFT h = rhs_complex.spec_div_element_wise(lhs_complex);

			DataArray<2> h_re = h.getRealWithDataArray();
			DataArray<2> h_im = h.getImagWithDataArray();

			DataArray<2> u0_dx_h_re = (io_u - op.diff_c_x(h_re));
			DataArray<2> u0_dx_h_im = op.diff_c_x(h_im);

			DataArray<2> v0_dy_h_re = (io_v - op.diff_c_y(h_re));
			DataArray<2> v0_dy_h_im = op.diff_c_y(h_im);

			Complex2DArrayFFT u0_dx_h(io_u.resolution);
			Complex2DArrayFFT v0_dy_h(io_u.resolution);

			u0_dx_h.loadRealAndImagFromDataArrays(u0_dx_h_re, u0_dx_h_im);
			v0_dy_h.loadRealAndImagFromDataArrays(v0_dy_h_re, v0_dy_h_im);

			Complex2DArrayFFT u = -alpha/kappa*u0_dx_h + i_parameters.sim_f/kappa*v0_dy_h;
			Complex2DArrayFFT v = -i_parameters.sim_f/kappa*u0_dx_h + -alpha/kappa*v0_dy_h;

			DataArray<2> u_re = u.getRealWithDataArray();
			DataArray<2> u_im = u.getImagWithDataArray();
			DataArray<2> v_re = v.getRealWithDataArray();
			DataArray<2> v_im = v.getImagWithDataArray();

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
