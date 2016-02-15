/*
 * RexiSWE_HelmholtzSolver.hpp
 *
 *  Created on: 2 Oct 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_PROGRAMS_REXISWE_REXISWE_HELMHOLTZSOLVER_HPP_
#define SRC_PROGRAMS_REXISWE_REXISWE_HELMHOLTZSOLVER_HPP_

#include <sweet/Complex2DArrayFFT.hpp>


class RexiSWE_HelmholtzSolver
{

	typedef
	bool (*smoother_fun)(
			std::complex<double> i_kappa,
			double i_gh0,
			const Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double *i_domain_size,
			double i_error_threshold,
			int i_max_iters,
			double i_omega,
			int i_verbosity
	);

public:


	/**
	 * Solve complex-valued Helmholtz problem
	 *
	 * (kappa - gh*D2) X = B
	 */
	static
	bool smoother_jacobi(
			std::complex<double> i_kappa,
			double i_gh0,
			const Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double *i_domain_size,
			double i_error_threshold = 0.000001,
			int i_max_iters = 999999999,
			double i_omega = 0.1,
			int i_verbosity = 0
	)
	{
#if 0
		std::cout << i_kappa << std::endl;
		std::cout << i_gh0 << std::endl;
//		std::cout << i_rhs << std::endl;
#endif

		double inv_helm_h[2];
		inv_helm_h[0] = (double)i_rhs.resolution[0]/(double)i_domain_size[0];
		inv_helm_h[1] = (double)i_rhs.resolution[1]/(double)i_domain_size[1];

		double scalar_Dx = inv_helm_h[0]*inv_helm_h[0];
		double scalar_Dy = inv_helm_h[1]*inv_helm_h[1];
		double scalar_C = -(2.0*(inv_helm_h[0]*inv_helm_h[0]) + 2.0*(inv_helm_h[1]*inv_helm_h[1]));

		std::cout << "kappa: " << i_kappa << std::endl;
		std::cout << "scalar_Dx: " << scalar_Dx << std::endl;
		std::cout << "scalar_Dy: " << scalar_Dy << std::endl;
		std::cout << "scalar_C: " << scalar_C << std::endl;

		std::complex<double> inv_diag = 1.0/(i_kappa - i_gh0*scalar_C);

//		double inv_diag_re = 1.0/(i_kappa.real() - i_gh0*scalar_C);
//		double inv_diag_im = 1.0/(i_kappa.imag());

		if (i_verbosity > 3)
			std::cout << "inv_diag: " << inv_diag << std::endl;

		Complex2DArrayFFT new_x(i_rhs.resolution);
		int i = 0;
		double residual = 0;
		for (i = 0; i < i_max_iters; i++)
		{
			Complex2DArrayFFT tmp = io_x;
#if 0
			Complex2DArrayFFT y = (i_rhs + i_gh0*io_x.op_stencil_Re_X(scalar_Dx, scalar_Dy));

			Complex2DArrayFFT y_re = y;
			y_re.setAllIm(0);

			Complex2DArrayFFT y_im = y;
			y_im.setAllRe(0);

			Complex2DArrayFFT t = y_re*inv_diag_re + y_im*inv_diag_im;
//			Complex2DArrayFFT t = (y_re+y_im)*inv_diag;
			io_x = (1.0-i_omega)*io_x + i_omega*t;

#else
			io_x = (1.0-i_omega)*io_x + i_omega*(
					(i_rhs + i_gh0*io_x.op_stencil_Re_X(scalar_Dx, scalar_Dy))*inv_diag
				);
#endif

			residual = ((io_x*i_kappa - i_gh0*io_x.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C)) - i_rhs).reduce_rms();

			if (i_verbosity > 3)
				std::cout << i << " RESIDUAL: " << residual << std::endl;

			if (residual < i_error_threshold)
				return true;
		}

		return false;
	}



	/**
	 * Solve complex-valued Helmholtz problem with CG solver
	 * It's not really a smoother, but it's nice to use the same interfaces as for the smoothers ;-).
	 *
	 * (kappa - gh*D2) X = B
	 *
	 * See "Numerical Analysis", page 474
	 */
	static
	bool smoother_conjugate_gradient(
			std::complex<double> i_kappa,
			double i_gh0,
			const Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double *i_domain_size,
			double i_error_threshold = 0.000001,
			int i_max_iters = 999999999,
			double i_omega = -1,		///< omega is non-sence for CG
			int i_verbosity = 0
	)
	{
		const Complex2DArrayFFT &b = i_rhs;
		Complex2DArrayFFT &x = io_x;

		double inv_helm_h[2];
		inv_helm_h[0] = (double)i_rhs.resolution[0]/(double)i_domain_size[0];
		inv_helm_h[1] = (double)i_rhs.resolution[1]/(double)i_domain_size[1];

		double scalar_Dx = 1.0*(inv_helm_h[0]*inv_helm_h[0]);
		double scalar_Dy = 1.0*(inv_helm_h[1]*inv_helm_h[1]);
		double scalar_C = -(2.0*(inv_helm_h[0]*inv_helm_h[0]) + 2.0*(inv_helm_h[1]*inv_helm_h[1]));

		if (	std::abs(scalar_C+i_kappa.real()) < std::abs(i_kappa.imag())	)
		{
			std::cerr << "ERROR: Diagonal dominance for CG solver not assured (TODO: check if this is required for CG) - stopping!" << std::endl;
			exit(-1);
		}

#define A(x)	(x*i_kappa - x.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C))

#if 0
		double foo = 1.0/std::sqrt(i_kappa.real() - i_gh0*scalar_C);
		double c_re = foo;
		double c_im = foo;

		c_re = 1;
		c_im = 1;

#	define Ci(x)	(x.multiply_real_imag(c_re, c_im))
#else
#	define Ci(x)	(x)
#endif

//#define Ci(x)	(x)

		Complex2DArrayFFT r = b - A(x);
		Complex2DArrayFFT w = Ci(r);
		Complex2DArrayFFT v = Ci(w);

		// TODO: should we maybe use the conjugate dot product?
		std::complex<double> alpha = w.dotProd(w);

		if (i_verbosity > 3)
			std::cout << "RESIDUAL: " << r.reduce_rms() << std::endl;


		int i = 0;
		for (i = 0; i < i_max_iters; i++)
		{
//			if (v.reduce_rms() < i_error_threshold)
//				return true;

			Complex2DArrayFFT u = A(v);
			std::complex<double> t = alpha / (v.dotProd(u));

			x = x + t*v;
			r = r - t*u;
			w = Ci(r);
			std::complex<double> beta = w.dotProd(w);

			if (std::abs(beta) < i_error_threshold)
			{
				if (r.reduce_rms() < i_error_threshold)
				{
					if (i_verbosity > 0)
						std::cout << "FIN RESIDUAL: " << (b-A(x)).reduce_rms() << " after " << i << " iterations" << std::endl;
					return true;
				}
			}

			if (i_verbosity > 3)
				std::cout << "RESIDUAL: " << r.reduce_rms() << std::endl;

			std::complex<double> s = beta / alpha;
			v = Ci(w) + v*s;

			alpha = beta;
		}

#undef A
#undef Ci

		return false;
	}



#if 0
	/**
	 * Solve complex-valued Helmholtz problem with CG solver
	 * It's not really a smoother, but it's nice to use the same interfaces as for the smoothers ;-).
	 *
	 * (kappa - gh*D2) X = B
	 *
	 * See "Numerical Analysis", page 474
	 */
	static
	bool smoother_conjugate_gradient_real(
			std::complex<double> i_kappa,
			double i_gh0,
			const Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			Operators2D &op,
			double *i_domain_size,
			double i_error_threshold = 0.000001,
			int i_max_iters = 999999999,
			double i_omega = -1,
			int i_verbosity = 0
	)
	{
		DataArray<2> b_re = i_rhs.toDataArrays_Real();
		DataArray<2> b_im = i_rhs.toDataArrays_Imag();

		DataArray<2> x_re = io_x.toDataArrays_Real();
		DataArray<2> x_im = io_x.toDataArrays_Imag();

#define A_re(x_re, x_im)	((x_re)*i_kappa.real() - (x_im)*(i_kappa.imag()) - i_gh0*(op.diff2_c_x(x_re) + op.diff2_c_y(x_re)))
#define A_im(x_re, x_im)	((x_re)*i_kappa.imag() + (x_im)*(i_kappa.real()) - i_gh0*(op.diff2_c_x(x_im) + op.diff2_c_y(x_im)))

		DataArray<2> r_re = b_re - A_re(x_re, x_im);
		DataArray<2> r_im = b_im - A_im(x_re, x_im);

		DataArray<2> v_re = r_re;
		DataArray<2> v_im = r_im;
		double alpha = (v_re*v_re).reduce_sum() + (v_im*v_im).reduce_sum();

		int i = 0;
		for (i = 0; i < i_max_iters; i++)
		{
			if (i_verbosity > 3)
			{
				DataArray<2> r_re = b_re - A_re(x_re, x_im);
				DataArray<2> r_im = b_im - A_im(x_re, x_im);

				std::cout << "RESIDUAL: " << std::sqrt(((r_re*r_re).reduce_sum() + (r_im*r_im).reduce_sum())/(double)(x_re.resolution[0]*x_re.resolution[1])) << std::endl;
			}

//			if (v.reduce_rms() < i_error_threshold)
//				return true;

			DataArray<2> u_re = A_re(v_re, v_im);
			DataArray<2> u_im = A_im(v_re, v_im);
			double t = alpha / ((v_re*u_re).reduce_sum() + (v_im*u_im).reduce_sum());

			x_re = x_re + t*v_re;
			x_im = x_im + t*v_im;

			r_re = r_re - t*u_re;
			r_im = r_im - t*u_im;

			double beta = (r_re*r_re).reduce_sum() + (r_im*r_im).reduce_sum();

			if (std::abs(beta) < i_error_threshold)
			{
				double res = std::sqrt(beta/(double)(x_re.resolution[0]*x_re.resolution[1]));
				if (res < i_error_threshold)
				{
					if (i_verbosity > 0)
						std::cout << "FIN RESIDUAL: " << res << " after " << i << " iterations" << std::endl;

					io_x.loadRealAndImagFromDataArrays(x_re, x_im);
					return true;
				}
			}

			double s = beta / alpha;
			v_re = r_re + v_re*s;
			v_im = r_im + v_im*s;

			alpha = beta;
		}

#undef A

		return false;
	}
#endif



	static
	bool multigrid(
			std::complex<double> i_kappa,
			double i_gh0,
			const Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			smoother_fun i_smoother_fun,
			double *i_domain_size,
			double i_error_threshold = 0.000001,
			int i_max_iters = 999999999,
			double i_omega = 0.1,
			int level_restriction = 2,
			int i_verbosity = 0
	)
	{
		int levels = 1;
		while ((std::size_t)(1 << levels) < std::min(i_rhs.resolution[0], i_rhs.resolution[1]))
			levels++;

		if (level_restriction < 0)
		{
			level_restriction = levels + level_restriction-1;
		}
		else if (level_restriction < 2)
		{
			std::cerr << "SAFETY STOP: going to 2x2 resolution generates instable solver, please restrict level" << std::endl;
			exit(1);
		}

		std::vector<Complex2DArrayFFT> rhs;
		std::vector<Complex2DArrayFFT> x;
		rhs.resize(levels-level_restriction);
		x.resize(levels-level_restriction);

		rhs[0] = i_rhs;
		x[0] = io_x;

		for (int l = 0; l < levels-1-level_restriction; l++)
		{
#if 0
			{
				std::stringstream ss;
				ss << "mg_" << l << "_x_re_a_down_pre.csv";
				x[l].toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l << "_x_im_a_down_pre.csv";
				x[l].toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}

			Complex2DArrayFFT residual = helmholtz_iterative_get_residual(i_kappa, i_gh0, rhs[l], x[l], i_domain_size);
			{
				std::stringstream ss;
				ss << "mg_" << l << "_residual_re_a_down_pre.csv";
				residual.toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l << "_residual_im_a_down_pre.csv";
				residual.toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}

			std::cout << "PRESMOOTH " << l << " with resolution " << rhs[l].resolution[0] << "x" << rhs[l].resolution[1] << std::endl;
			// presmoother
			std::cout << "           RESIDUAL START " << helmholtz_iterative_get_residual_rms(i_kappa, i_gh0, rhs[l], x[l], i_domain_size) << std::endl;
#endif

#if 1
			i_smoother_fun(
					i_kappa,
					i_gh0,
					rhs[l],		// rhs
					x[l],		// x
					i_domain_size,
					i_error_threshold,
					1,
					i_omega,
					i_verbosity
			);
#endif
			if (i_verbosity > 2)
				std::cout << "           RESIDUAL END " << helmholtz_iterative_get_residual_rms(i_kappa, i_gh0, rhs[l], x[l], i_domain_size) << std::endl;

#if 0
			{
				std::stringstream ss;
				ss << "mg_" << l << "_x_re_b_down_post.csv";
				x[l].toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l << "_x_im_b_down_post.csv";
				x[l].toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}

			residual = helmholtz_iterative_get_residual(i_kappa, i_gh0, rhs[l], x[l], i_domain_size);
			{
				std::stringstream ss;
				ss << "mg_" << l << "_residual_re_b_down_post.csv";
				residual.toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l << "_residual_im_b_down_post.csv";
				residual.toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}
#endif

			rhs[l+1] = rhs[l].scale_down();
			x[l+1] = x[l].scale_down();
		}

		{
			int l = levels-1-level_restriction;
			assert(l > 0);

			if (i_verbosity > 2)
				std::cout << "MID SMOOTH " << l << " with resolution " << rhs[l].resolution[0] << "x" << rhs[l].resolution[1] << std::endl;

			bool retval = i_smoother_fun(
					i_kappa,
					i_gh0,
					rhs[l],
					x[l],
					i_domain_size,
					i_error_threshold,
					i_max_iters,
					i_omega,
					i_verbosity
			);

			if (!retval)
			{
				std::cout << "No convergence for iterative solver" << std::endl;
				return false;
			}
		}

		for (int l = levels-1-level_restriction; l >= 1; l--)
		{
			if (i_verbosity > 2)
				std::cout << "           RESIDUAL PREV LEVEL CHECK " << helmholtz_iterative_get_residual_rms(i_kappa, i_gh0, rhs[l], x[l], i_domain_size) << std::endl;

			x[l-1] = x[l].scale_up();

#if 0
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_x_re_c_up_pre.csv";
				x[l-1].toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_x_im_c_up_pre.csv";
				x[l-1].toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}

			Complex2DArrayFFT residual = helmholtz_iterative_get_residual(i_kappa, i_gh0, rhs[l-1], x[l-1], i_domain_size);
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_residual_re_c_up_pre.csv";
				residual.toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_residual_im_c_up_pre.csv";
				residual.toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}
#endif

#if 0
			std::cout << std::cout.precision(4);
			std::cout << x[l] << std::endl;
			std::cout << rhs[l] << std::endl;
			std::cout << std::endl;
			std::cout << x[l-1] << std::endl;
			std::cout << rhs[l-1] << std::endl;
			exit(-1);
#endif

			if (i_verbosity > 2)
			{
				std::cout << "POSTSMOOTH " << l-1 << " with resolution " << rhs[l-1].resolution[0] << "x" << rhs[l-1].resolution[1] << std::endl;
				std::cout << "           RESIDUAL START " << helmholtz_iterative_get_residual_rms(i_kappa, i_gh0, rhs[l-1], x[l-1], i_domain_size) << std::endl;
			}

			bool retval = i_smoother_fun(
					i_kappa,
					i_gh0,
					rhs[l-1],	// rhs
					x[l-1],		// x
					i_domain_size,
					i_error_threshold,
					i_max_iters,
					i_omega,
					i_verbosity
			);

			if (i_verbosity > 2)
				std::cout << "           RESIDUAL END " << helmholtz_iterative_get_residual_rms(i_kappa, i_gh0, rhs[l-1], x[l-1], i_domain_size) << std::endl;

#if 0
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_x_re_d_up_post.csv";
				x[l-1].toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_x_im_d_up_post.csv";
				x[l-1].toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}


			residual = helmholtz_iterative_get_residual(i_kappa, i_gh0, rhs[l-1], x[l-1], i_domain_size);
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_residual_re_d_up_post.csv";
				residual.toDataArrays_Real().file_saveData_ascii(ss.str().c_str());
			}
			{
				std::stringstream ss;
				ss << "mg_" << l-1 << "_residual_im_d_up_post.csv";
				residual.toDataArrays_Imag().file_saveData_ascii(ss.str().c_str());
			}
#endif

			if (!retval)
			{
				std::cout << "No convergence for iterative solver" << std::endl;
				return false;
			}
		}

		io_x = x[0];
		return true;
	}


	static
	double helmholtz_iterative_get_residual_rms(
			std::complex<double> i_kappa,
			double i_gh0,
			Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double *i_domain_size
	)
	{
		double helm_h[2] =
					{
						(double)i_domain_size[0] / (double)i_rhs.resolution[0],
						(double)i_domain_size[1] / (double)i_rhs.resolution[1]
					};

		std::complex<double> D2_kernel[3*3];
		{
			D2_kernel[0] = { 0, 0 };
			D2_kernel[1] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[2] = { 0, 0 };

			D2_kernel[3] = { 1.0/(helm_h[0]*helm_h[0]), 0 };
			D2_kernel[4] = { -(2.0/(helm_h[0]*helm_h[0]) + 2.0/(helm_h[1]*helm_h[1])), 0};
			D2_kernel[5] = { 1.0/(helm_h[0]*helm_h[0]), 0 };

			D2_kernel[6] = { 0, 0 };
			D2_kernel[7] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[8] = { 0, 0 };
		}

		return ((io_x*i_kappa - i_gh0*io_x.op_stencil_3x3(D2_kernel)) - i_rhs).reduce_rms();
	}



	static
	Complex2DArrayFFT helmholtz_iterative_get_residual(
			std::complex<double> i_kappa,
			double i_gh0,
			Complex2DArrayFFT &i_rhs,
			Complex2DArrayFFT &io_x,
			double *i_domain_size
	)
	{
		double helm_h[2] =
					{
						(double)i_domain_size[0] / (double)i_rhs.resolution[0],
						(double)i_domain_size[1] / (double)i_rhs.resolution[1]
					};

		std::complex<double> D2_kernel[3*3];
		{
			D2_kernel[0] = { 0, 0 };
			D2_kernel[1] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[2] = { 0, 0 };

			D2_kernel[3] = { 1.0/(helm_h[0]*helm_h[0]), 0 };
			D2_kernel[4] = { -(2.0/(helm_h[0]*helm_h[0]) + 2.0/(helm_h[1]*helm_h[1])), 0};
			D2_kernel[5] = { 1.0/(helm_h[0]*helm_h[0]), 0 };

			D2_kernel[6] = { 0, 0 };
			D2_kernel[7] = { 1.0/(helm_h[1]*helm_h[1]), 0 };
			D2_kernel[8] = { 0, 0 };
		}

		return ((io_x*i_kappa - i_gh0*io_x.op_stencil_3x3(D2_kernel)) - i_rhs);
	}
};



#endif /* SRC_PROGRAMS_REXISWE_REXISWE_HELMHOLTZSOLVER_HPP_ */
