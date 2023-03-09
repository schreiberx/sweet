/*
 * PlaneOperatorsComplex.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_PLANE_OPERATORS_COMPLEX_HPP_
#define SRC_INCLUDE_SWEET_PLANE_OPERATORS_COMPLEX_HPP_


#include <sweet/core/plane/PlaneData_Config.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>

namespace sweet
{

class PlaneOperatorsComplex
{
public:
	ErrorBase error;
	
	PlaneData_Config *planeDataConfig;

	// differential operators
	PlaneData_SpectralComplex diff_c_x, diff_c_y;
	PlaneData_SpectralComplex diff2_c_x, diff2_c_y;

	/**
	 * D2, e.g. for viscosity
	 */
	inline PlaneData_SpectralComplex diff2(
			const PlaneData_SpectralComplex &i_dataArray
	)
	{
		return diff2_c_x(i_dataArray) + diff2_c_y(i_dataArray);
	}



	/**
	 *        __2
	 * apply  \/  operator (aka Laplace)
	 */
	inline PlaneData_SpectralComplex laplace(
			const PlaneData_SpectralComplex &i_a
	)
	{
		return diff2_c_x(i_a)+diff2_c_y(i_a);
	}


	/**
	 *        __
	 * apply  \/ .  operator
	 */
	inline PlaneData_SpectralComplex diff_dot(
			const PlaneData_SpectralComplex &i_a
	)
	{
		return diff_c_x(i_a)+diff_c_y(i_a);
	}



	/**
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	inline PlaneData_SpectralComplex diffN_x(
			const PlaneData_SpectralComplex &io_u,
			int i_order
	)
	{
		if (i_order == 0)
			return io_u;

		PlaneData_SpectralComplex tu = io_u;

		for (int i = 0; i < i_order/2; i++)
			tu = diff2_c_x(tu);

		if (i_order & 1)
			tu = diff_c_x(tu);

		return tu;
	}


	/**
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	inline PlaneData_SpectralComplex diffN_y(
			const PlaneData_SpectralComplex &io_v,
			int i_order
	)
	{
		if (i_order == 0)
			return io_v;

		PlaneData_SpectralComplex tv = io_v;

		for (int i = 0; i < i_order/2; i++)
			tv = diff2_c_y(tv);

		if (i_order & 1)
			tv = diff_c_y(tv);

		return tv;
	}


#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Diffusion or hyperviscosity coefficients
	 * Simply calculates the spectral coefficients
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 *
	 * Returns operator D^q
	 *
	 */
	inline PlaneData_SpectralComplex diffusion_coefficient(
			int i_order
	)
	{
		//Check if even
		assert( i_order % 2 == 0);
		assert( i_order > 0);
		PlaneData_SpectralComplex out = diff2_c_x+diff2_c_y;

		for (int i = 1; i < i_order/2; i++)
			out = pow(-1, i)*(diff2_c_x(out)+diff2_c_y(out));

		return out;
	}

	/**
	 * Calculates implicit diffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	inline PlaneData_SpectralComplex implicit_diffusion(
			const PlaneData_SpectralComplex &i_data,
			double i_coef,
			int i_order
	)
	{
		PlaneData_SpectralComplex out=i_data;


		PlaneData_SpectralComplex asdf = operator*(1.3, out);

		// Get diffusion coefficients (these are the -mu*dt*D^q, where q is the order
		PlaneData_SpectralComplex diff = -i_coef*diffusion_coefficient(i_order);

		// Add 1 to get denominator
		diff = diff.spectral_addScalarAll(1.0);

		// Invert
		diff = diff.spectral_invert();
		// apply to data
		out=diff(out);
		return out;
	}
#endif


	bool setup(
			PlaneData_Config *i_planeDataConfig,		///< data config setup for spectral transformations
			ShackPlaneDataOps *i_shackPlaneDataOps
	)
	{
		return setup(
				i_planeDataConfig,
				i_shackPlaneDataOps->plane_domain_size
		);
	}


	bool setup(
			PlaneData_Config &i_planeDataConfig,		///< data config setup for spectral transformations
			ShackPlaneDataOps &i_shackPlaneDataOps
	)
	{
		return setup(
				&i_planeDataConfig,
				i_shackPlaneDataOps.plane_domain_size
		);
	}

	bool setup(
			PlaneData_Config &i_planeDataConfig,		///< data config setup for spectral transformations
			ShackPlaneDataOps *i_shackPlaneDataOps
	)
	{
		return setup(
				&i_planeDataConfig,
				i_shackPlaneDataOps->plane_domain_size
		);
	}

	bool setup(
			PlaneData_Config *i_planeDataConfig,
			const double i_domain_size[2]
	)
	{
		assert(planeDataConfig == nullptr);
		planeDataConfig = i_planeDataConfig;


		/*
		 * Setup spectral differential operators
		 * 		diff(e(ix), x)
		 */
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
		std::cerr << "Activate spectral space during compile time to use spectral diffs. Otherwise, the convolution would be freakingly expensive" << std::endl;
		assert(false);
#endif

		/*
		 * DIFF X
		 */
		{
			diff_c_x.setup(planeDataConfig);
			diff_c_x.spectral_set_zero();
			double scale_x = 2.0*M_PI/i_domain_size[0];

			/*
			 * left bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[0][1][0]; j < planeDataConfig->spectral_complex_ranges[0][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[0][0][0]; i < planeDataConfig->spectral_complex_ranges[0][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, (double)i*scale_x);

			/*
			 * left top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[1][1][0]; j < planeDataConfig->spectral_complex_ranges[1][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[1][0][0]; i < planeDataConfig->spectral_complex_ranges[1][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, (double)i*scale_x);

			/*
			 * right bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[2][1][0]; j < planeDataConfig->spectral_complex_ranges[2][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[2][0][0]; i < planeDataConfig->spectral_complex_ranges[2][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, -(double)(planeDataConfig->spectral_complex_ranges[2][0][1]-i)*scale_x);

			/*
			 * right top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[3][1][0]; j < planeDataConfig->spectral_complex_ranges[3][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[3][0][0]; i < planeDataConfig->spectral_complex_ranges[3][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, -(double)(planeDataConfig->spectral_complex_ranges[3][0][1]-i)*scale_x);
		}

		/*
		 * DIFF Y
		 */
		{
			diff_c_y.setup(planeDataConfig);

			diff_c_y.spectral_set_zero();
			double scale_y = 2.0*M_PI/i_domain_size[1];


			/*
			 * left bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[0][1][0]; j < planeDataConfig->spectral_complex_ranges[0][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[0][0][0]; i < planeDataConfig->spectral_complex_ranges[0][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, (double)j*scale_y);


			/*
			 * left top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[1][1][0]; j < planeDataConfig->spectral_complex_ranges[1][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[1][0][0]; i < planeDataConfig->spectral_complex_ranges[1][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, -(double)(planeDataConfig->spectral_complex_ranges[1][1][1]-j)*scale_y);

			/*
			 * right bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[2][1][0]; j < planeDataConfig->spectral_complex_ranges[2][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[2][0][0]; i < planeDataConfig->spectral_complex_ranges[2][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, (double)(j)*scale_y);

			/*
			 * right top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = planeDataConfig->spectral_complex_ranges[3][1][0]; j < planeDataConfig->spectral_complex_ranges[3][1][1]; j++)
				for (std::size_t i = planeDataConfig->spectral_complex_ranges[3][0][0]; i < planeDataConfig->spectral_complex_ranges[3][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, -(double)(planeDataConfig->spectral_complex_ranges[3][1][1]-j)*scale_y);
		}

		/*
		 * 2nd order differential operators
		 */
		/*
		 * DIFF2 X
		 */
		{
			diff2_c_x.setup(planeDataConfig);

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < planeDataConfig->spectral_complex_array_data_number_of_elements; i++)
				diff2_c_x.spectral_space_data[i] = diff_c_x.spectral_space_data[i]*diff_c_x.spectral_space_data[i];
		}


		/*
		 * DIFF2 X
		 */
		{
			diff2_c_y.setup(planeDataConfig);

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < planeDataConfig->spectral_complex_array_data_number_of_elements; i++)
				diff2_c_y.spectral_space_data[i] = diff_c_y.spectral_space_data[i]*diff_c_y.spectral_space_data[i];
		}

		return true;
	}



public:
	PlaneOperatorsComplex()	:
		planeDataConfig(nullptr)
	{
	}



public:
	PlaneOperatorsComplex(
		PlaneData_Config *i_planeDataConfig,
		const double i_domain_size[2]	///< domain size
	)	:
		planeDataConfig(i_planeDataConfig),

		diff_c_x(i_planeDataConfig),
		diff_c_y(i_planeDataConfig),
		diff2_c_x(i_planeDataConfig),
		diff2_c_y(i_planeDataConfig)
	{
		setup(i_planeDataConfig, i_domain_size);
	}
};

}

#endif
