/*
 * GSLMultiMin.hpp
 *
 * Compute the extrema of an arbitrary dimensional function
 *
 *  Created on: Jan 30, 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_LIBMATH_GSLMULTIMIN_HPP_
#define SRC_INCLUDE_LIBMATH_GSLMULTIMIN_HPP_


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multimin.h>


template <int D, typename T>
class GSLMultiMin
{
	const gsl_multimin_fminimizer_type *minimizer_type;
	gsl_multimin_fminimizer *minimizer;
	gsl_multimin_function gminex_func;

	gsl_vector *iter_coord;
	gsl_vector *iter_step_size;


	/**
	 * coordinate of local maximum extrema
	 *
	 * the components of this array are updated according to the current 1D quadrature
	 */
public:
	std::array<T,D> coord;

	typedef T (*fun_template)(
			const T *i_coords,
			const void *i_param
	);

	/**
	 * function to run quadrature for
	 */
	fun_template fun;

	/**
	 * user parameters for function
	 */
	const void *fun_user_param;

	/**
	 * start and end of integration interval
	 */
	std::array<T,D> start, end;


public:
	GSLMultiMin(
			int i_quad_key = 6
	)	:
		fun(nullptr),
		fun_user_param(nullptr)
	{
		minimizer_type = gsl_multimin_fminimizer_nmsimplex2;
		minimizer = gsl_multimin_fminimizer_alloc(minimizer_type, D);

		iter_coord = gsl_vector_alloc(D);
		iter_step_size = gsl_vector_alloc(D);
	}



public:
	~GSLMultiMin()
	{
		gsl_multimin_fminimizer_free(minimizer);

		gsl_vector_free(iter_coord);
		gsl_vector_free(iter_step_size);
	}



	/**
	 * This class wraps the execution of the function
	 * run_quadrature_dim
	 * into a class to allow executing the quadrature in parallel
	 */
	class RunQuadratureObjectWrapper
	{
	public:
		GSLMultiMin *gslMultiMin;

		static
		double
		run_multimin(double *x, std::size_t dim, void *i_param)
		{
			AssertE(dim == D);
			RunQuadratureObjectWrapper *r = (RunQuadratureObjectWrapper*)i_param;
			return r->gslQuadrature->fun(x, r->gslQuadrature->fun_user_param);
		}
	};



	static
	double funCallback(
			const gsl_vector *v,
			void *params
	)
	{
		GSLMultiMin &m = *(GSLMultiMin*)params;
		T coords[D];

		// setup coordinates and clamp to area
		for (int i = 0; i < D; i++)
		{
			coords[i] = gsl_vector_get(v, i);

			if ((coords[i] < m.start[i]) || (coords[i] > m.end[i]))
				return std::numeric_limits<T>::max();

//			coords[i] = std::min(std::max(gsl_vector_get(v, i), m.start[i]), m.end[i]);
		}

		return m.fun(coords, m.fun_user_param);
	}


	/**
	 * run the multiminimizer
	 *
	 * \return true on success
	 */
private:
	T run_multimin()
	{
		gminex_func.n = D;
		gminex_func.f = funCallback;
		gminex_func.params = this;

		gsl_vector_set_all(iter_step_size, 0.001);

		gsl_multimin_fminimizer_set(
				minimizer,
				&gminex_func,
				iter_coord,			///< starting point
				iter_step_size		///< step size
			);

#if 0
		std::cout << std::endl;
		std::cout << "ITER: " << std::endl;
#endif
		int status;
		std::size_t iter = 0;
		do
		{
#if 0
			std::cout << iter << ": ";
			for (int i = 0; i < D; i++)
				std::cout << gsl_vector_get(minimizer->x, i) << ", ";
			std::cout << std::endl;
#endif
			status = gsl_multimin_fminimizer_iterate(minimizer);

			if (status)
			{
				std::cout << "FAIL" << std::endl;
				return -std::numeric_limits<T>::max();
			}

			double size = gsl_multimin_fminimizer_size(minimizer);
//			std::cout << "size: " << size << std::endl;
			status = gsl_multimin_test_size(size, (end[0]-start[0])*1e-3);

			iter++;
		}
		while (status == GSL_CONTINUE && iter < 10000);

//		AssertE(status == GSL_SUCCESS);
#if 0
		std::cout << "SUCCESS (" << iter << ") ";
		for (int i = 0; i < D; i++)
			coord[i] = gsl_vector_get(minimizer->x, i);
		std::cout << "COORD: " << coord << std::endl;
#endif

		return -funCallback(minimizer->x, this);
	}


	/**
	 * run quadrature over function with the interval given by its center and size
	 */
public:
	T multimin_center_size(
			fun_template i_fun,			///< function to compute quadrature for
			const void *i_fun_user_param,		///< user-defined parameters
			const std::array<T,D> &i_center,	///< center of interval
			const std::array<T,D> &i_size		///< size of interval
	)
	{
		fun = i_fun;
		fun_user_param = i_fun_user_param;

		// initialize interval
		for (int i = 0; i < D; i++)
		{
			start[i] = i_center[i] - i_size[i]*(T)0.5;
			end[i] = i_center[i] + i_size[i]*(T)0.5;
			gsl_vector_set(iter_coord, i, i_center[i]);
		}

		return run_multimin();
	}



	/**
	 * run quadrature over function with the interval given by lower (start) and higher (end) coordinate
	 */
public:
	T multimin_start_end(
			fun_template i_fun,			///< function to compute quadrature for
			const void *i_fun_user_param,		///< user-defined parameters
			const std::array<T,D> &i_start,	///< start of interval
			const std::array<T,D> &i_end		///< end of interval
	)
	{
		fun = i_fun;
		fun_user_param = i_fun_user_param;

		start = i_start;
		end = i_end;

		for (int i = 0; i < D; i++)
			gsl_vector_set(iter_coord, i, (i_start[i]+i_end[i])*(T)0.5);

		return run_multimin();
	}
};



#endif /* SRC_INCLUDE_LIBMATH_GSLQUADRATURE_HPP_ */
