/*
 * RexiSWE_HelmholtzSolver.hpp
 *
 *  Created on: 2 Oct 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_PROGRAMS_REXI_REXISWE_HELMHOLTZSOLVER_HPP_
#define SRC_PROGRAMS_REXI_REXISWE_HELMHOLTZSOLVER_HPP_

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

//		std::cout << "scalar_Dx: " << scalar_Dx << std::endl;
//		std::cout << "scalar_Dy: " << scalar_Dy << std::endl;
//		std::cout << "scalar_C: " << scalar_C << std::endl;

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
#if 0
			if (i == 1778)
			{
				Complex2DArrayFFT asdf = io_x-tmp;

				std::cout << std::setprecision(4);
				std::cout << tmp << std::endl;
				std::cout << std::endl;
				std::cout << io_x << std::endl;
				std::cout << std::endl;
				std::cout << asdf << std::endl;
				std::cout << std::endl;

				residual = ((tmp*i_kappa - i_gh0*tmp.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C)) - i_rhs).reduce_rms();
				std::cout << "RESIDUAL pre: " << residual << std::endl;

				residual = ((io_x*i_kappa - i_gh0*io_x.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C)) - i_rhs).reduce_rms();
				std::cout << "RESIDUAL: " << residual << std::endl;

				std::cout << (io_x*i_kappa - i_gh0*io_x.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C) - i_rhs) << std::endl;

				Complex2DArrayFFT k = (io_x*i_kappa - i_gh0*io_x.op_stencil_Re_X_C(scalar_Dx, scalar_Dy, scalar_C) - i_rhs);

				double t = 0;
				for (int j = 0; j < io_x.resolution[1]; j++)
					for (int i = 0; i < io_x.resolution[0]; i++)
						t += k.getRe(j, i)*k.getRe(j, i);

				double l = 0;
				for (int j = 0; j < io_x.resolution[1]; j++)
					for (int i = 0; i < io_x.resolution[0]; i++)
						l += k.getIm(j, i)*k.getIm(j, i);

				double u = 0;
				for (int j = 0; j < io_x.resolution[1]; j++)
					for (int i = 0; i < io_x.resolution[0]; i++)
						u += k.getRe(j, i)*k.getRe(j, i) + k.getIm(j, i)*k.getIm(j, i);

				std::cout << t << std::endl;
				std::cout << l << std::endl;
				std::cout << u << std::endl;

				double sum = std::sqrt(u/double(k.resolution[0]*k.resolution[1]));
				std::cout << sum << std::endl;
				std::cout << k.reduce_rms_quad() << std::endl;

				for (int j = 0; j < io_x.resolution[1]; j++)
					for (int i = 0; i < io_x.resolution[0]; i++)
						if (std::abs(asdf.get(j, i)) > 1e-7)
							std::cout << j << ", " << i << ", " << asdf.get(j, i) << std::endl;

				std::cout << std::endl;
				exit(-1);
			}
#endif
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
		double alpha = (w*w).reduce_sum_real_imag();

		int i = 0;
		for (i = 0; i < i_max_iters; i++)
		{
			if (i_verbosity > 3)
				std::cout << "RESIDUAL: " << (b-A(x)).reduce_rms() << std::endl;

//			if (v.reduce_rms() < i_error_threshold)
//				return true;

			Complex2DArrayFFT u = A(v);
			double t = alpha / (v*u).reduce_sum_real_imag();

			x = x + t*v;
			r = r - t*u;
			w = Ci(r);
			double beta = (w*w).reduce_sum_real_imag();

			if (std::abs(beta) < i_error_threshold)
			{
				if (r.reduce_rms() < i_error_threshold)
				{
					if (i_verbosity > 0)
						std::cout << "FIN RESIDUAL: " << (b-A(x)).reduce_rms() << " after " << i << " iterations" << std::endl;
					return true;
				}
			}

			double s = beta / alpha;
			v = Ci(w) + v*s;

			alpha = beta;
		}

#undef A
#undef Ci

		return false;
	}



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



#endif /* SRC_PROGRAMS_REXI_REXISWE_HELMHOLTZSOLVER_HPP_ */
