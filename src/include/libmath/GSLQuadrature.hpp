/*
 * GSLQuadrature.hpp
 *
 *  Created on: Jan 30, 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_LIBMATH_GSLQUADRATURE_HPP_
#define SRC_INCLUDE_LIBMATH_GSLQUADRATURE_HPP_


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>


#define GSL_RANDOMIZED_QUADRATURE	1


template <int D, typename T>
class GSLQuadrature
{
#if GSL_RANDOMIZED_QUADRATURE
	gsl_monte_plain_state *mc_s;
	const gsl_rng_type *mc_trng;
	gsl_rng *mc_rng;

#else

	/**
	 * size of workspace
	 */
	const std::size_t workspace_size = 1024;


	gsl_integration_workspace *workspace;


	/**
	 * coordinate for high-dimensional integration
	 *
	 * the components of this array are updated according to the current 1D quadrature
	 */
	std::array<T,D> coord;

#endif

	typedef T (*fun_template)(const T *i_coords, const void *i_param);

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

	// error thresholds for gsl quadrature function
	T epsabs = 0;
	T epsrel = 1e-8;
	int quad_key;


public:
	GSLQuadrature(
			int i_quad_key = 6
	)	:
		fun(nullptr),
		fun_user_param(nullptr),
		quad_key(i_quad_key)
	{
#if GSL_RANDOMIZED_QUADRATURE
	    mc_s = gsl_monte_plain_alloc(D);
		gsl_rng_env_setup();
		mc_trng = gsl_rng_default;
		mc_rng = gsl_rng_alloc(mc_trng);
#else
		workspace = gsl_integration_workspace_alloc(workspace_size);
#endif
	}



public:
	~GSLQuadrature()
	{
#if GSL_RANDOMIZED_QUADRATURE
	    gsl_monte_plain_free(mc_s);
	    gsl_rng_free(mc_rng);
#else
		gsl_integration_workspace_free(workspace);
#endif
	}



	/**
	 * This class wraps the execution of the function
	 * run_quadrature_dim
	 * into a class to allow executing the quadrature in parallel
	 */
	class RunQuadratureObjectWrapper
	{
	public:
		GSLQuadrature *gslQuadrature;

#if GSL_RANDOMIZED_QUADRATURE

		static
		double
		run_quadrature_dim(double *x, std::size_t dim, void *i_param)
		{
			AssertE(dim == D);
			RunQuadratureObjectWrapper *r = (RunQuadratureObjectWrapper*)i_param;
			return r->gslQuadrature->fun(x, r->gslQuadrature->fun_user_param);
		}

#else
		int dim;

		static
		double
		run_quadrature_dim(double x, void *i_param)
		{
			RunQuadratureObjectWrapper *r = (RunQuadratureObjectWrapper*)i_param;
			return r->gslQuadrature->run_quadrature_dim(x, r->dim);
		}
#endif
	};



#if !GSL_RANDOMIZED_QUADRATURE
	/**
	 * run the quadrature
	 */
private:
	double run_quadrature_dim(
		T x,			///< quadrature point
		int i_dim		///< dimension of quadrature point
	)
	{
		coord[i_dim] = x;
		if (i_dim == (D-1))
		{
			// last dimension: evaluate function
			return fun(coord.data(), fun_user_param);
		}

		// next dimension to run the quadrature for
		RunQuadratureObjectWrapper r;
		r.gslQuadrature = this;
		r.dim = i_dim+1;

		double result, error;

		gsl_function F;
		F.function = &RunQuadratureObjectWrapper::run_quadrature_dim;
		F.params = &r;

		// use version which includes singularities due to abs() function
		gsl_integration_qag(
				&F,					///< function
				start[r.dim], end[r.dim],	///< interval
				epsabs,				///< epsabs
				epsrel,				///< epsrel
				workspace_size,		///< limit
				quad_key,
				workspace,
				&result,			///< result
				&error				///< absolute error
			);

		return result;
	}
#endif


private:
	T run_quadrature()
	{
#if GSL_RANDOMIZED_QUADRATURE
		// next dimension to run the quadrature for
		RunQuadratureObjectWrapper r;
		r.gslQuadrature = this;

		gsl_monte_function F;
		F.f = &RunQuadratureObjectWrapper::run_quadrature_dim;
		F.dim = D;
		F.params = &r;

		int calls = std::max(500, (int)((end[0]-start[0])*500000.0));
//		std::cout << calls << std::endl;

		double result, error;
	    gsl_monte_plain_integrate(
	    		&F,
				start.data(),
				end.data(),
				D,
	    		calls,
				mc_rng,
				mc_s,
				&result,
				&error
			);

	    return result;

#else
		epsabs = 1e-7;
		epsrel = 1e-9;

		// next dimension to run the quadrature for
		RunQuadratureObjectWrapper r;
		r.gslQuadrature = this;
		r.dim = 0;

		gsl_function F;
		F.function = &RunQuadratureObjectWrapper::run_quadrature_dim;
		F.params = &r;

		double result, error;

		// use version which includes singularities due to abs() function
		gsl_integration_qag(
				&F,					///< function
				start[0], end[0],	///< interval
				epsabs,				///< epsabs
				epsrel,				///< epsrel
				workspace_size,		///< limit
				quad_key,
				workspace,
				&result,			///< result
				&error				///< absolute error
			);

		return result;
#endif

	}


	/**
	 * run quadrature over function with the interval given by its center and size
	 */
public:
	T quadrature_center_size(
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
		}

		return run_quadrature();
	}



	/**
	 * run quadrature over function with the interval given by lower (start) and higher (end) coordinate
	 */
public:
	T quadrature_start_end(
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

		return run_quadrature();
	}


};



#endif /* SRC_INCLUDE_LIBMATH_GSLQUADRATURE_HPP_ */
