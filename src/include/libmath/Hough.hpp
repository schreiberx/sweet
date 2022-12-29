/*
 * Hough.hpp
 *
 *  Created on: 20 Dec 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_HOUGH_HPP_
#define SRC_INCLUDE_HOUGH_HPP_

#include <cmath>
#include <limits>
#include <iostream>
#include <functional>
#include <sweet/SWEETError.hpp>
#include <Eigen/Eigenvalues>


class ALPCoefficients
{
public:
	int size;
	double* A = nullptr;
	double* B = nullptr;

	ALPCoefficients(int i_max_degree)
	{
		this->size = Hough::sizeP(i_max_degree);
		this->A = MemBlockAlloc::alloc<double>(size * sizeof(double));
		this->B = MemBlockAlloc::alloc<double>(size * sizeof(double));
	}

}

class Hough
{
public:

	double g  = 9.8;
	double h0 = 1e3;
	double a = 2e7 / M_PI;
	double omega = 2. * M_PI / (60. * 60. * 24.);

	double gamma = std::sqrt(g * h0) / ( 2. * a * omega);
	double epsl  = 1. / (gamma * gamma);

	ALPCoefficients* coeffs = nullptr;

	Hough(int i_max_degree)
	{
		this->coeffs = new ALPCoefficients(i_max_degree);
	}

	//sizeP(maxDegree)
	//Return the size of the set of Associated Legendre Polynomials ``P_l^m(x)`` of
	//degree less than or equal to the given maximum degree
	int sizeP(int i_max_degree)
	{
		return int( (i_max_degree + 1) * (i_max_degree + 2) / 2);
	}

	//sizeY(maxDegree)
	//Return the size of the set of real spherical harmonics ``Y_{l,m}(θ,φ)`` of·
	//degree less than or equal to the given maximum degree
	int sizeY(int o_max_degree)
	{
		return (i_max_degree + 1) * (i_max_degree + 1);
	}

	//index_p(l,m)
	//Return the index into a flat array of Associated Legendre Polynomials ``P_l^m``
	//for the given indices ``(l,m)``.
	//``P_l^m`` are stored in l-major order i.e. [P(0,0), [P(1,0), P(1,1), P(2,0), ...]
	int index_p(int i_l, int i_m)
	{
		return i_m + int(i_l * (i_l + 1) / 2) + 1
	}

	//index_y(l,m)
	//Return the index into a flat array of real spherical harmonics ``Y_{l,m}``
	//for the given indices ``(l,m)``.
	//``Y_{l,m}`` are stored in l-major order i.e.·
	//[Y(0,0), [Y(1,-1), Y(1,0), Y(1,1), Y(2,-2), ...]
	int index_y(int i_l, int i_m)
	{
		return i_m + i_l + (i_l * i_l) + 1
	}






	double* legendreP(int i_l, int i_m, double i_x, int i_lmax)
	{
		int ind = this->index_p(i_l, i_m);
		return std::sqrt(M_PI) * this->compute_p(i_lmax, i_x)[ind];
	}


	double* legendreP(int i_l, int i_m, double* i_x, int i_lmax, int size)
	{
		double* y = nullptr;
		y = MemBlockAlloc::alloc<double>(size * sizeof(double));
		for (int i = 0; i < size; i++)
			y[i] = this->legendreP(i_l, i_m, i_x[i], i_lmax)
		return y;
	}


	double* legendreP(int i_lmax, double* i_x)
	{
		return std::sqrt(M_PI) * this->compute_p(i_lmax, i_x);
	}

	ALPCoefficients* compute_coefficients(int i_L)
	{
		for (int l = 2; l <= L; l++)
		{
			int ls = l * l;
			int lm1s = (l - 1) * (l - 1);
			for (int m = 0; m <= l - 2; m++)
			{
				int ms = m * m;
				int idx = this->index_p(l, m);
				this->coeffs->A[idx] = std::sqrt( (4. * ls - 1.) / (ls - ms) );
				this->coeffs->B[idx] = -std::sqrt( (lm1s - ms) / (4. * lm1s - 1.));
			}
		}
		return this->coeffs;
	}



/////!·······compute_coefficients(L)
/////Create an array large enough to store an entire set of Associated Legendre·
/////Polynomials ``P_l^m(x)`` of maximum degree L.·
/////"""
/////function allocate_p(L::Int)
/////!·······P = Array{Float64}(undef,sizeP(L))
/////end


//compute_p(L, x, coeff, P)
//Compute an entire set of Associated Legendre Polynomials ``P_l^m(x)``
//using the given coefficients, and store in the array P.
	void compute_P(int i_L, double i_x, ALPCoefficients* i_coeffs, double* io_P)
	{
		int size = this->sizeP(L);
		assert(i_coeffs.size >= size);

		double sintheta = std::sqrt(1. - i_x * i_x);
		double temp = istd::sqrt(.5 / M_PI);
		P[this->index_p(0, 0)] = temp;

		if (i_L > 0)
		{
			P[this->index_p(1, 0)] = i_x * std::sqrt(3) * temp;
			temp = -std::sqrt(1.5) * sintheta * temp;
			P[this->index_p(1, 1)] = temp;

			for (int l = 2; l <= L; l++)
			{
				for (int m = 0; m <= l-2; m++)
				{
					int idx = this->index_p(l, m) =   i_coeffs.A[idx] * (i_x * io_P[this->index_p(l - 1, m)]
									+ i_coeffs.B[idx] * io_P[this->idxed_p(l - 2, m)] );
				}
				io_P[this->index_p(l, l - 1)] = i_x * std::sqrt(2. * (l - 1) + 3) * temp;
				temp = -std::sqrt(1. + .5 / l) * sintheta * temp;
				io_P[this->index_p(l, l)] = temp;
			}
		}

	}


//compute_p(L, x)
//Compute an entire set of Associated Legendre Polynomials ``P_l^m(x)`` where
//``0 ≤ l ≤ L`` and ``0 ≤ m ≤ l``. Assumes ``|x| ≤ 1``.
	double* compute_p(int i_L, double i_x)
	{
		double* P = nullptr;
		int size = this->sizeP(i_L);
		P = MemBlockAlloc::alloc<double>(size * sizeof(double));
		ALPCoefficients* coeff = new ALPCoefficients(i_L);
		this->compute_P(i_L, i_x, coeff, P);
		return P;
	}


//compute_y(L, P, φ, Y)
//Compute an entire set of real spherical harmonics ``Y_{l,m}(θ,φ)``·
//using the given Associated Legendre Polynomials ``P_l^m(cos θ)``
//and store in array Y
	void compute_y(int i_L, double* i_P, double i_phi, double* o_Y)
	{
		for (int l = 0; l <= i_L; l++)
			o_Y[this->index_y(l, 0)] = P[this->index_p(l, 0)] * .5 * std::sqrt(2);
	}

	double c1 = 1.;
	double c2 = std::cos(i_phi);
	double s1 = 0.;
	double s2 = -std::sin(i_phi);
	double tc = 2. * c2;

	for (int m = 1; m <= L; m++)
	{
		double s = tc * s1 - s2;
		double c = tc * c1 - c2;
		s2 = s1;
		s1 = s;
		c2 = c1;
		c1 = c;
		for (int l = m; l <= L; l++)
		{
			o_Y[this->index_y(l, -m)] = P[this->index_p(l, m)] * s;
			o_Y[this->index_y(l, m)] = P[this->index_p(l, m)] * c;
		}
	}


	//compute_y(L, x, φ)
	//Compute an entire set of real spherical harmonics ``Y_{l,m}(θ,φ)`` for·
	//``x = cos θ`` where ``0 ≤ l ≤ L`` and ``-l ≤ m ≤ l``.
	double* compute_y(int i_L, double i_x, double i_phi)
	{
		int size = this->sizeP(i_L);
		P = MemBlockAlloc::alloc<double>(size * sizeof(double));
		ALPCoefficients* coeff = new ALPCoefficients(i_L);
		this->compute_p(i_L, i_x, coeff, P);

		size = this->sizeY(L);
		Y = MemBlockAlloc::alloc<double>(size * sizeof(double));
		this->compute_y(i_L, P, i_phi, Y);
		return Y;
	}


//////////	double* compute_p(int i_lmax, double* i_x)
//////////	{
//////////		double* pmn = nullptr;
//////////		int size = 
//////////		pmn = MemBlockAlloc::alloc<double>(size * sizeof(double));
//////////	}
//////////function compute_p(lmax::Int64, x::Array{Float64,1})
//////////        pmn=Array{Float64}(undef, length(x), sizeP(lmax))
//////////        for (xi, xval) in enumerate(x)
//////////                pmn[xi, :] = compute_p(lmax, xval)
//////////        end
//////////        return pmn
//////////end


	double e2(int i_n, int i_m)
	{
		double aux1 = i_n * i_n - i_m * i_m;
		double aux2 = 4. * i_n * i_n - 1;
		return std::sqrt(aux1 / aux2);
	}

	double q2(int i_n, int i_m)
	{
		return i_n / (i_n + 1) * this->e2(i_n + 1, i_m);
	}

	double p2(int i_n, int i_m)
	{
		return (i_n + 1) / i_n * this->e2(i_n, i_m);
	}

	double k2(int i_n, int i_m)
	{
		return -i_m / (i_n * (i_n + 1.));
	}

	double r(int i_n)
	{
		return -i_n * (i_n + 1) * gamma * gamma;
	}

	Eigen::MatrixXd pmn_matrix_a(int i_s, int i_n)
	{
		int nn = 3 * i_n;
		int lm = int(nn / 3 - 1);
		Eigen::MatrixXd a(nn, nn);
		for (int i = 1; i <= lm; i++)
		{
			int m = int(3 * i - 1);
			double t1 = i_s + 2. * (i - 1.);
			double t2 = t1 + 1.;
			a(m - 1,m-1 - 1)   = this->q2(t2-1,i_s);
			a(m - 1,m - 1)     = this->k2(t2,i_s);
			a(m - 1,m+2 - 1)   = this->p2(t2+1,i_s);
			a(m+1 - 1,m-1 - 1) = this->r(t2-1);
			a(m+2 - 1,m - 1)   = this->q2(t2,i_s);
			a(m+2 - 1,m+2 - 1) = this->k2(t2+1,i_s);
			a(m+2 - 1,m+3 - 1) = this->p2(t2+2,i_s);
			a(m+2 - 1,m+4 - 1) = -1;
		}
		a(1 - 1,1 - 1) = this->k2(i_s,i_s);
		a(1 - 1,2 - 1) = this->p2(s+1,i_s);
		a(1 - 1,3 - 1) = -1.;
		double t1 = i_s + 2. * (nn / 3. -1.);
		double t2 = t1 + 1;
		a(nn-1 - 1,nn-2 - 1) = this->q2(t1,i_s);
		a(nn - 1,nn-2 - 1) = this->r(t1);
		a(nn-1 - 1,nn-1 - 1) = this->k2(t2,i_s);
		return a
	}

	Eigen::MatrixXd pmn_matrix_b(int i_s, int i_n)
	{
		int nn = 3 * i_n;
		int lm = int(nn / 3 - 1);
		Eigen::MatrixXd a(nn, nn);
		for (int i = 1; i <= lm; i++)
		{
			int m = int(3 * i - 1);
			double t1 = i_s + 2. * (i - 1.);
			double t2 = t1 + 1.;
			a(m - 1,m-1 - 1)   = this->q2(t1,i_s);
			a(m - 1,m - 1)     = this->k2(t2,i_s);
			a(m - 1,m+1 - 1)   = -1.;
			a(m - 1,m+2 - 1)   = this->p2(t2+1,i_s);
			a(m+1 - 1,m - 1)   = this->r(t2);
			a(m+2 - 1,m - 1)   = this->q2(t2,i_s);
			a(m+2 - 1,m+2 - 1) = this->k2(t2+1,i_s);
			a(m+2 - 1,m+3 - 1) = this->p2(t2+2,i_s);
		}
		a(1 - 1,1 - 1) = this->k2(i_s,i_s);
		a(1 - 1,2 - 1) = this->p2(i_s+1,i_s);
		double t1 = s + 2. * (nn / 3. -1.);
		double t2 = t1 + 1;
		a(nn-1 - 1,nn-2 - 1) = this->q2(t1,i_s);
		a(nn-1 - 1,nn-1 - 1) = this->k2(t2,i_s);
		a(nn-1 - 1,nn - 1) = -1;
		a(nn - 1,nn-1 - 1) = this->r(t1);
		return a
	}


	double freq(int i_s, int i_l, int i_n, int i_kind)
	{
		if (i_kind == 2)
		{
			Eigen::MatrixXd b;
			if (i_ll%2 == 1)
			{
				b = this->pmn_matrix_a(i_s, i_n);
				i_l = int((i_l-1)/2 + 1);
			}
			else
			{
				b = this->pmn_matrix_b(i_s, i_n);
				i_l = int((i_l)/2 + 1);
			}
			Eigen::VectorXcd eig = b.eigenvalues();
			return eig(i_n + i_l);
		}
		else
		{
			Eigen::MatrixXd b;
			if (i_ll%2 == 0)
			{
				b = this->pmn_matrix_a(i_s, i_n);
				i_l = int((i_l)/2 + 1);
			}
			else
			{
				b = this->pmn_matrix_b(i_s, i_n);
				i_l = int((i_l - 1)/2 + 1);
			}
			Eigen::VectorXcd eig = b.eigenvalues();
			if (i_kind == 1)
				return eig(i_n - 1 + i_l); // to check
			else if (i_kind == 3)
				return eig(2 * i_n + i_l); // to check
		}
	}


	bool check_sym(int i_s. int i_l, int i_n, int i_kind)
	{
		if (i_kind == 2)
		{
			if (l%2 == 1)
				return true;
			else
				return false;
		}
		else
		{
			if (l%2 == 0)
				return true;
			else
				return false;
		}
	}

	Eigen::VectorXd eigen_vec(int i_s, int i_l, int i_n, int i_kind)
	{
		bool sym = this->check_sym(i_s, i_l, i_n, i_kind);
		double freq_foo = this->freq(i_s, i_l, i_n, i_kind);
		Eigen::MatrixXd matrix;
		if (sym)
			matrix = this->pmn_matrix_a(i_s, i_n);
		else
			matrix = this->pmn_matrix_a(i_s, i_n);
		Eigen::EigenSolver<Eigen::MatrixXd> eig_solver(matrix);

		if (i_l % 2 == 0)
			i_l = int(i_l / 2 + 1);
		else
			i_l = int((i_l - 1) / 2 + 1);
		int idx;
		if (kind == 1)
			idx = i_n + 1 - i_l;
		else if (kind == 2)
			idx = i_n + i_l;
		else if (kind == 3)
			idx = 2 * i_n + i_l;
		/// supposing eigenvalues are sorted!!!
		/////ind = sortperm(eig_val)[ind_foo]
		return eig_solver.eigenvectors().col(idx);
	}



!·······return u, -1im.*v, z,  du, -1im.*dv, dz,  1im*s.*u, s.*v, 1im*s.*z, dive, vrt

	void uvz_and_deriv(	int i_s, int i_l, int i_lmax, int i_kind, double i_x, int i_nlat,
				double* o_u,
				double* o_miv,
				double* o_z,
				double* o_du,
				double* o_midv,
				double* o_dz,
				double* o_pisu,
				double* o_sv,
				double* o_pisz,
				double* o_div,
				double* o_vort
			)
	{
		Eigen::VectorXd foo = this->eigen_vec(i_s, i_l, i_lmax, i_kind);
		bool sym = this->check_sym(i_s, i_l, i_lmax, i_kind);
		int lmax2 = i_lmax * i_lmax;
		int lmax2s = lmax2 + i_s;

		double pmn = this->legendreP(lmax2s, i_x);

		double* coef_a = nullptr;
		double* coef_b = nullptr;
		double* coef_c = nullptr;

		coef_a = MemBlockAlloc::alloc<double>(lmax2 * sizeof(double));
		coef_b = MemBlockAlloc::alloc<double>(lmax2 * sizeof(double));
		coef_c = MemBlockAlloc::alloc<double>(lmax2 * sizeof(double));

		// Gathering/ordering expansion coefficients
		for (int i = 1; i <= i_lmax; i++)
		{
			if (sym)
			{
				coef_a[2*i-1 - 1] = foo[3*i-2 - 1];
				coef_b[2*i - 1]   = foo[3*i-1 - 1];
				coef_c[2*i-1 - 1] = foo[3*i - 1];
			}
			else
			{
				coef_b[2*i-1 - 1] = foo[3*i-2 - 1];
				coef_a[2*i - 1]   = foo[3*i-1 - 1];
				coef_c[2*i - 1]   = foo[3*i - 1];
			}
		}

		// Calculating normalizing constant
		double const_K = 0;
		for (int i = 1; i <= lmax2; i++)
		{
			int m = i_s + i - 1;
			const_K += (coef_a[i - 1] * coef_a[i - 1] + coef_b[i - 1] * coef_b[i - 1]) * (m * (m + 1)) / epsl + coef_c[i - 1] * coef_c[i - 1];
		}
		const_K = std::sqrt(const_K);
		// Normalizing coefficients
		for (int i = 0; i < lmax2; i++)
		{
			coef_a[i] /= const_K;
			coef_b[i] /= const_K;
			coef_c[i] /= const_K;
		}

		// Calculating psi, phi, z, dpsi, dphi
		std::vector<std::vector<std::complex<double>>> xi = {};
		for (int i = 0; i < 7; i++)
		{
			std::vector<std::complex<double>> aux(i_nlat, 1i);
			xi.push_back(aux);
		}

		for (int i = 1; i <= lmax2; i++)
		{
			int m = i_s + i - 1;
			lpmn = pmn[:, this->index_p(i_m, i_s)];
			for (int j = 0; j < nlat; j++)
			{
				xi[0, j] += 1i * coef_a[i - 1] * lpmn[j];
				xi[1, j] += coef_b[i - 1] * lpmn[j];
				xi[2, j] += coef_c[i - 1] * lpmn[j];
			}

			double* lpmn_m1 = nullptr;
			lpmn_m1 = MemBlockAlloc::alloc<double>(nlat * sizeof(double));
			if (m == 1)
				for (int j = 0; j < nlat; j++)
					lpmn_m1[j] = 0.;
			else
			{
				int idx = this->index_p(m - 1, i_s);
				for (int j = 0; j < nlat; j++)
					lpmn_m1[j] = pmn[j, idx];
			}

			double* lpmn_p1 = nullptr;
			lpmn_p1 = MemBlockAlloc::alloc<double>(nlat * sizeof(double));
			if (m == lmax2)
				for (int j = 0; j < nlat; j++)
					lpmn_p1[j] = 0.;
			else
			{
				int idx = this->index_p(m + 1, i_s);
				for (int j = 0; j < nlat; j++)
					lpmn_p1[j] = pmn[j, idx];
			}

			double* dlpmn = nullptr;
			dlpmn = MemBlockAlloc::alloc<double>(nlat * sizeof(double));
			double e21 = this->e2(m, i_s);
			double e22 = this->e2(m + 1, i_s);
			for (int j = 0; j < nlat; j++)
				dplmn[j] = (m + 1) * e21 * lpmn_m1[j] - m * e22 * lpmn_p1[j];

			for (int j = 0; j < nlat; j++)
			{
				xi[3, j] += coef_b[i - 1] * dplmn[j];
				xi[4, j] += 1i * coef_a[i - 1] * dplmn[j];
				xi[5, j] += -m * (m + 1) * 1i * coef_a[i - 1] * lpmn[j];
				xi[6, j] += -m * (m + 1) * coef_b[i - 1] * lpmn[j];
			}

		}

		// Calculating u, v, z
		for (int j = 0; j < nlat; j++)
		{
			u[j] = (1i * i_s * xi[0, j] - xi[3, j]) / std::sqrt(epsl * 1. - i_x * i_x);
			v[j] = 1i * (1i * i_s * xi[1, j] + xi[4, j]) / std::sqrt(epsl * 1. - i_x * i_x);
			z[j] = xi[2, j];

			div[j] = xi[5, j] / std::sqrt(epsl);
			vrt[j] = xi[6, j] / std::sqrt(epsl);
		}

		// Calculating du, dv, dz
		for (int j = 0; j < nlat; j++)
		{
			du[j] = i_s * v[j];
			dv[j] = i_s * u[j];
			dz[j] = 0.;
		}

		for (int i = 1; i <= lmax2; i++)
		{
			int m = i_s + i - 1;
			lpmn = pmn[:, this->index_p(i_m, i_s)];

			// Checking position on lpmn values matrix
			double* lpmn_m1 = nullptr;
			lpmn_m1 = MemBlockAlloc::alloc<double>(nlat * sizeof(double));
			if (m == 1)
				for (int j = 0; j < i_nlat; j++)
					lpmn_m1[j] = 0.;
			else
			{
				int idx = this->index_p(m - 1, i_s);
				for (int j = 0; j < i_nlat; j++)
					lpmn_m1[j] = pmn[j, idx];
			}

			double* lpmn_p1 = nullptr;
			lpmn_p1 = MemBlockAlloc::alloc<double>(nlat * sizeof(double));
			if (m == lmax2)
				for (int j = 0; j < i_nlat; j++)
					lpmn_p1[j] = 0.;
			else
			{
				int idx = this->index_p(m + 1, i_s);
				for (int j = 0; j < i_nlat; j++)
					lpmn_p1[j] = pmn[j, idx];
			}

			for (int j = 0; j < i_nlat; j++)
			{
				du[j] += std::sqrt( (1. - i_x * i_x) / epsl) * m * (m + 1.) * coef_b[i - 1] * lpmn[j];
				dv[j] += std::sqrt( (1. - i_x * i_x) / epsl) * m * (m + 1.) * coef_a[i - 1] * lpmn[j];
				dz[j] += 1. / std::sqrt(1. - i_x * i_x) * coef_c[i - 1] * ( (m + 1) * this->e2(m, s) * lpmn_m1[j] - m * this->e2(m + 1, s) * lpmn_p1[j]);
			}

		}


		for (int j = 0; j < i_nlat; j++)
		{
			o_u[j] = u[j];
			o_miv[j] = -1i * v[j];
			o_z[j] = z[j];
			o_du[j] = du[j];
			o_midv[j] = -1i * dv[j];
			o_pisu[j] = 1i * i_s * u[j];
			o_sv[j] = i_s * v[j];
			o_pisz[j] = 1i * i_s * z[j];
			o_div = div[j];
			o_vrt = vrt[j];
		}


	}

};

