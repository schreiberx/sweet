/*
 * EigenMatrix3.hpp
 *
 *  Created on: 21 Nov 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *      Calculates the eigenvalues of a real or std::complex 3x3 matrices
 *
 *   ***** DO NOT USE******** CONFLICTS IN TYPES *********
 *
 */

#ifndef SRC_INCLUDE_EIGENMATRIX3_HPP_
#define SRC_INCLUDE_EIGENMATRIX3_HPP_

#include <cmath>
#include <limits>
#include <iostream>
#include <functional>
#include <libmath/DQStuff.hpp>
#include <sweet/SWEETError.hpp>
#include <complex>

//using namespace std;

class EigenMatrix3
{
public:
public:
	double matrix[3][3];	///< storage for data


	 /******************* CONSTRUCTORS *********************/

	/**
	 * constructor: load identity matrix
	 */
	inline EigenMatrix3()
	{
		SWEETError("DO NOT USE EigenMatrix3.hpp : CONFLICTS IN TYPES IN FINAL FORMULA ");
		loadIdentity();
	}


	/**
	 * constructor: load matrix entries
	 */
	inline EigenMatrix3(
			double v00,
			double v01,
			double v02,
			double v10,
			double v11,
			double v12,
			double v20,
			double v21,
			double v22
	)
	{
		SWEETError("DO NOT USE EigenMatrix3.hpp : CONFLICTS IN TYPES IN FINAL FORMULA ");
		matrix[0][0] = v00;
		matrix[0][1] = v01;
		matrix[0][2] = v02;

		matrix[1][0] = v10;
		matrix[1][1] = v11;
		matrix[1][2] = v12;

		matrix[2][0] = v20;
		matrix[2][1] = v21;
		matrix[2][2] = v22;

	}

	/*Functions*/

	/**
	 * load identity matrix
	 */
	inline EigenMatrix3& loadIdentity()
	{
		SWEETError("DO NOT USE EigenMatrix3.hpp : CONFLICTS IN TYPES IN FINAL FORMULA ");

		matrix[0][0] = (double)1;
		matrix[0][1] = (double)0;
		matrix[0][2] = (double)0;

		matrix[1][0] = (double)0;
		matrix[1][1] = (double)1;
		matrix[1][2] = (double)0;

		matrix[2][0] = (double)0;
		matrix[2][1] = (double)0;
		matrix[2][2] = (double)1;

		return *this;
	}


	std::complex<double> l_sqrtcplx(const std::complex<double> &i_value)
				{
		return std::sqrt(i_value);
				};

	std::complex<double> l_cbrtcplx(const std::complex<double> &i_value)
				{
		return std::pow(i_value,1.0/3.0);
				};

	void eigen3real(
			std::complex<double> o_eval[3]
	)
	{
		SWEETError("DO NOT USE EigenMatrix3.hpp : CONFLICTS IN TYPES IN FINAL FORMULA ");

#if 0
		/*
		double a00 = matrix[0][0];
		double a01 = matrix[0][1];
		double a02 = matrix[0][2];
		double a10 = matrix[1][0];
		double a11 = matrix[1][1];
		double a12 = matrix[1][2];
		double a20 = matrix[2][0];
		double a21 = matrix[2][1];
		double a22 = matrix[2][2];
		*/

		std::complex<double> I;
		I = -1;
		I = l_sqrtcplx(std::complex<double>(I));

		//Obtained via sympy using
		//   sweet/doc/misc/eigen/complexeigen3x3.py

		//o_eval[0] = (1.0L/6.0L)*((2*a00 + 2*a11 + 2*a22 - std::pow(2, 2.0L/3.0L)* l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 + l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) \
				+ std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) \
				- 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)))*l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 \
				+ l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 \
				- 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)) \
				+ 2*l_cbrtcplx(2)*(3*a00*a11 + 3*a00*a22 - 3*a01*a10 - 3*a02*a20 + 3*a11*a22 - 3*a12*a21 - std::pow(a00 + a11 + a22, 2)))/l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 \
				+ l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 \
				- 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21));
	    //o_eval[1] = (1.0L/12.0L)*((1 - l_sqrtcplx(3)*I)*(4*a00 + 4*a11 + 4*a22 + std::pow(2, 2.0L/3.0L)*(1 - l_sqrtcplx(3)*I)*l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 + l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)))*l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 + l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)) + 8*l_cbrtcplx(2)*(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2)))/((1 - l_sqrtcplx(3)*I)*l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 + l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)));
		//o_eval[2] = (1.0L/12.0L)*((1 + l_sqrtcplx(3)*I)*(4*a00 + 4*a11 + 4*a22 + std::pow(2, 2.0L/3.0L)*(1 + l_sqrtcplx(3)*I)*l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 + l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)))*l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 + l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)) + 8*l_cbrtcplx(2)*(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2)))/((1 + l_sqrtcplx(3)*I)*l_cbrtcplx(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 + l_sqrtcplx(-4*std::pow(-3*a00*a11 - 3*a00*a22 + 3*a01*a10 + 3*a02*a20 - 3*a11*a22 + 3*a12*a21 + std::pow(a00 + a11 + a22, 2), 3) + std::pow(-27*a00*a11*a22 + 27*a00*a12*a21 + 27*a01*a10*a22 - 27*a01*a12*a20 - 27*a02*a10*a21 + 27*a02*a11*a20 - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21), 2)) - 2*std::pow(a00 + a11 + a22, 3) + 9*(a00 + a11 + a22)*(a00*a11 + a00*a22 - a01*a10 - a02*a20 + a11*a22 - a12*a21)));
#endif

		o_eval[0]=1;
		o_eval[1]=1;
		o_eval[2]=1;

		return;
	}




};

#endif /* SRC_INCLUDE_EIGENMATRIX3_HPP_ */
