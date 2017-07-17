/*
 * PlaneOperators.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_INCLUDE_SWEET_PLANEOPERATORS_HPP_
#define SRC_INCLUDE_SWEET_PLANEOPERATORS_HPP_


#if SWEET_USE_PLANE_SPECTRAL_SPACE
//	#include <sweet/plane/PlaneDataComplex.hpp>
#endif

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataConfig.hpp>


class PlaneOperators
{
public:
	PlaneDataConfig *planeDataConfig;

public:
	// differential operators (central / forward / backward)
	PlaneData diff_c_x, diff_c_y;
	PlaneData diff_f_x, diff_f_y;
	PlaneData diff_b_x, diff_b_y;

	PlaneData diff2_c_x, diff2_c_y;

	PlaneData avg_f_x, avg_f_y;
	PlaneData avg_b_x, avg_b_y;

	PlaneData shift_left;
	PlaneData shift_right;
	PlaneData shift_up;
	PlaneData shift_down;


	/**
	 * D2, e.g. for viscosity
	 */
	PlaneData diff2(
			const PlaneData &i_dataArray
	)
	{
		return diff2_c_x(i_dataArray) + diff2_c_y(i_dataArray);
	}



	/**
	 *        __2
	 * apply  \/  operator (aka Laplace)
	 */
	inline PlaneData laplace(
			const PlaneData &i_a
	)
	{
		return diff2_c_x(i_a)+diff2_c_y(i_a);
	}


	/**
	 *        __
	 * apply  \/ .  operator
	 */
	inline PlaneData diff_dot(
			const PlaneData &i_a
	)
	{
		return diff_c_x(i_a)+diff_c_y(i_a);
	}



	/**
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	inline PlaneData diffN_x(
			const PlaneData &io_u,
			int i_order
	)
	{
		if (i_order == 0)
			return io_u;

		PlaneData tu = io_u;

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
	inline PlaneData diffN_y(
			const PlaneData &io_v,
			int i_order
	)
	{
		if (i_order == 0)
			return io_v;

		PlaneData tv = io_v;

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
	inline PlaneData diffusion_coefficient(
			int i_order
	)
	{
		//Check if even
		assert( i_order % 2 == 0);
		assert( i_order > 0);
		PlaneData out = diff2_c_x+diff2_c_y;

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
	inline PlaneData implicit_diffusion(
			const PlaneData &i_data,
			double i_coef,
			int i_order
	)
	{
		PlaneData out=i_data;

		// Get diffusion coefficients (these are the -mu*dt*D^q, where q is the order
		PlaneData diff = -i_coef*diffusion_coefficient(i_order);

		// Add 1 to get denominator
		diff = diff.spectral_addScalarAll(1.0);

		// Invert
		diff = diff.spectral_invert();
		// apply to data
		out=diff(out);
		return out;
	}
#endif

	PlaneOperators()	:
		planeDataConfig(nullptr),

		diff_c_x(1),
		diff_c_y(1),

		diff_f_x(1),
		diff_f_y(1),
		diff_b_x(1),
		diff_b_y(1),

		diff2_c_x(1),
		diff2_c_y(1),

		avg_f_x(1),
		avg_f_y(1),
		avg_b_x(1),
		avg_b_y(1),

		shift_left(1),
		shift_right(1),
		shift_up(1),
		shift_down(1)
	{

	}


	void setup(
		PlaneDataConfig *i_planeDataConfig,		///< data config setup for spectral transformations
		const double i_domain_size[2],			///< domain size
		bool i_use_spectral_basis_diffs = true	///< use spectral differentiation (d/dx e^ix)
	)
	{
		planeDataConfig = i_planeDataConfig;

		diff_c_x.setup(i_planeDataConfig);
		diff_c_y.setup(i_planeDataConfig);

		diff_f_x.setup(i_planeDataConfig);
		diff_f_y.setup(i_planeDataConfig);
		diff_b_x.setup(i_planeDataConfig);
		diff_b_y.setup(i_planeDataConfig);

		diff2_c_x.setup(i_planeDataConfig);
		diff2_c_y.setup(i_planeDataConfig);

		avg_f_x.setup(i_planeDataConfig);
		avg_f_y.setup(i_planeDataConfig);
		avg_b_x.setup(i_planeDataConfig);
		avg_b_y.setup(i_planeDataConfig);

		shift_left.setup(i_planeDataConfig);
		shift_right.setup(i_planeDataConfig);
		shift_up.setup(i_planeDataConfig);
		shift_down.setup(i_planeDataConfig);

		setup(i_domain_size, i_use_spectral_basis_diffs);
	}


	void setup(
			const double i_domain_size[2],
			bool i_use_spectral_basis_diffs
	)
	{

		double h[2] = {
				(double)i_domain_size[0] / (double)planeDataConfig->physical_res[0],
				(double)i_domain_size[1] / (double)planeDataConfig->physical_res[1]
		};

/////////////////////////////////////////////////////////////////////

		double avg_f_x_kernel[3][3] = {
				{0,0,0},
				{0,1,1},
				{0,0,0},
		};
		avg_f_x.kernel_stencil_setup(avg_f_x_kernel, 0.5);

		double avg_f_y_kernel[3][3] = {
				{0,1,0},
				{0,1,0},
				{0,0,0},
		};
		avg_f_y.kernel_stencil_setup(avg_f_y_kernel, 0.5);

		double avg_b_x_kernel[3][3] = {
				{0,0,0},
				{1,1,0},
				{0,0,0},
		};
		avg_b_x.kernel_stencil_setup(avg_b_x_kernel, 0.5);

		double avg_b_y_kernel[3][3] = {
				{0,0,0},
				{0,1,0},
				{0,1,0},
		};
		avg_b_y.kernel_stencil_setup(avg_b_y_kernel, 0.5);

/////////////////////////////////////////////////////////////////////

		double shift_left_kernel[3][3] = {
				{0,0,0},
				{0,0,1},
				{0,0,0},
		};
		shift_left.kernel_stencil_setup(shift_left_kernel);

		double shift_right_kernel[3][3] = {
				{0,0,0},
				{1,0,0},
				{0,0,0},
		};
		shift_right.kernel_stencil_setup(shift_right_kernel);

		double shift_up_kernel[3][3] = {
				{0,0,0},
				{0,0,0},
				{0,1,0},
		};
		shift_up.kernel_stencil_setup(shift_up_kernel);

		double shift_down_kernel[3][3] = {
				{0,1,0},
				{0,0,0},
				{0,0,0},
		};
		shift_down.kernel_stencil_setup(shift_down_kernel);

/////////////////////////////////////////////////////////////////////

		if (i_use_spectral_basis_diffs)
		{
			/*
			 * setup spectral differential operators
			 * 		diff(e(ix), x)
			 */
			// Assume, that errors are linearly depending on the resolution
			// see test_spectral_ops.cpp

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
			std::cerr << "Activate spectral space during compile time to use spectral diffs. Otherwise, the convolution would be freakingly expensive" << std::endl;
			assert(false);
			exit(-1);
#else

			/*
			 * Note, that there's a last column which is set to 0 (Nyquist freq, noise in signal)
			 * PXT: removed this setting to zero (changed < to <=), because of 2nd and higher order differentiation
			 * MaS: changed <= to < for the x-axis because of invalid memory access...
			 */
			diff_c_x.spectral_set_all(0, 0);

			for (int j = 0; j <= (int)planeDataConfig->spectral_data_size[1]/2; j++)
			{
				for (int i = 0; i < (int)planeDataConfig->spectral_data_size[0]; i++)
				{
					std::size_t idxa = ((j)*planeDataConfig->spectral_data_size[0])+(i);
					diff_c_x.spectral_space_data[idxa] = {0.0, (double)(((double)i)*2.0*M_PIl/(double)i_domain_size[0])};

					std::size_t idxb = ((diff_c_x.planeDataConfig->spectral_data_size[1]-1-j)*planeDataConfig->spectral_data_size[0])+(i);
					diff_c_x.spectral_space_data[idxb] = {0.0, (double)((double)i*2.0*M_PIl/(double)i_domain_size[0])};

#if 0
					if (i == (int)planeDataConfig->spectral_data_size[0]-1)
					{
						diff_c_x.spectral_space_data[idxa] *= 2.0;
						diff_c_x.spectral_space_data[idxb] *= 2.0;
					}
#endif
				}
			}


			/*
			 * DIFF operator in y axis
			 */
			diff_c_y.spectral_set_all(0, 0);
			// TODO: shift j for loop by +1
			for (int j = 0; j <= (int)diff_c_y.planeDataConfig->spectral_data_size[1]/2-1; j++)
			{
				for (int i = 0; i < (int)diff_c_y.planeDataConfig->spectral_data_size[0]; i++)
				{
					std::size_t idxa = ((j+1)*planeDataConfig->spectral_data_size[0])+(i);
					diff_c_y.spectral_space_data[idxa] = {0, (double)((double)(j+1)*2.0*M_PIl/(double)i_domain_size[1])};

					std::size_t idxb = ((diff_c_y.planeDataConfig->spectral_data_size[1]-(j+1))*planeDataConfig->spectral_data_size[0])+(i);
					diff_c_y.spectral_space_data[idxb] = {0, (double)(-(double)(j+1)*2.0*M_PIl/(double)i_domain_size[1])};

#if 0
					diff_c_y.spectral_set(
							j+1, i,
							0,
							(double)(j+1)*2.0*M_PIl/(double)i_domain_size[1]
						);
					diff_c_y.spectral_set(
							diff_c_y.planeDataConfig->spectral_data_size[1]-(j+1), i,
							0,
							-(double)(j+1)*2.0*M_PIl/(double)i_domain_size[1]
						);
#endif
				}
			}


			/**
			 * TODO: WARNING! These operators are setup in Cartesian space,
			 * hence they are not as accurate as spectral operators
			 */
			double d_f_x_kernel[3][3] = {
					{0,0,0},
					{0,-1,1},
					{0,0,0}
			};
			diff_f_x.kernel_stencil_setup(d_f_x_kernel, 1.0/h[0]);

			double d_f_y_kernel[3][3] = {
					{0,1,0},
					{0,-1,0},
					{0,0,0},
			};
			diff_f_y.kernel_stencil_setup(d_f_y_kernel, 1.0/h[1]);


			double d_b_x_kernel[3][3] = {
					{0,0,0},
					{-1,1,0},
					{0,0,0}
			};
			diff_b_x.kernel_stencil_setup(d_b_x_kernel, 1.0/h[0]);

			double d_b_y_kernel[3][3] = {
					{0,0,0},
					{0,1,0},
					{0,-1,0},
			};
			diff_b_y.kernel_stencil_setup(d_b_y_kernel, 1.0/h[1]);



#if 1
			/*
			 * 2nd order differential operators
			 */
			diff2_c_x.spectral_set_all(0, 0);

			for (int j = 0; j <= (int)diff_c_x.planeDataConfig->spectral_data_size[1]/2; j++)
			{
				for (int i = 0; i < (int)diff_c_x.planeDataConfig->spectral_data_size[0]; i++)
				{
					double fac = (double)i*2.0*M_PIl/(double)i_domain_size[0];

					std::size_t idxa = ((j)*planeDataConfig->spectral_data_size[0])+(i);
					diff2_c_x.spectral_space_data[idxa] = {-fac*fac, 0.0};

					std::size_t idxb = ((diff_c_x.planeDataConfig->spectral_data_size[1]-1-j)*planeDataConfig->spectral_data_size[0])+(i);
					diff2_c_x.spectral_space_data[idxb] = {-fac*fac, 0.0};

#if 0
					diff2_c_x.spectral_set(
							j, i,
							-fac*fac,
							0
						);
					diff2_c_x.spectral_set(
							diff_c_x.planeDataConfig->spectral_data_size[1]-1-j, i,
							-fac*fac,
							0
						);
#endif
				}
			}


			diff2_c_y.spectral_set_all(0, 0);

			for (int j = 0; j <= (int)diff_c_y.planeDataConfig->spectral_data_size[1]/2-1; j++)
			{
				for (int i = 0; i < (int)diff_c_y.planeDataConfig->spectral_data_size[0]; i++)
				{
					double fac = (double)(j+1)*2.0*M_PIl/(double)i_domain_size[1];

					std::size_t idxa = ((j+1)*planeDataConfig->spectral_data_size[0])+(i);
					diff2_c_y.spectral_space_data[idxa] = {-fac*fac, 0.0};

					std::size_t idxb = ((diff_c_y.planeDataConfig->spectral_data_size[1]-(j+1))*planeDataConfig->spectral_data_size[0])+(i);
					diff2_c_y.spectral_space_data[idxb] = {-fac*fac, 0.0};

#if 0
					diff2_c_y.spectral_set(
							j+1, i,
							-fac*fac,
							0
						);

					diff2_c_y.spectral_set(
							diff_c_y.planeDataConfig->spectral_data_size[1]-(j+1), i,
							-fac*fac,
							0
						);
#endif
				}
			}
#else
			diff2_c_x = diff_c_x(diff_c_x);
			diff2_c_y = diff_c_y(diff_c_y);
#endif

#endif
		}
		else
		{
			double diff1_x_kernel[3][3] = {
					{0,0,0},
					{-1.0,0,1.0},
					{0,0,0}
			};
			diff_c_x.kernel_stencil_setup(diff1_x_kernel, 1.0/(2.0*h[0]));

			double diff1_y_kernel[3][3] = {
					{0,1.0,0},	// higher y coordinate
					{0,0,0},
					{0,-1.0,0},	// lower y coordinate
			};
			diff_c_y.kernel_stencil_setup(diff1_y_kernel, 1.0/(2.0*h[1]));

			double d_f_x_kernel[3][3] = {
					{0,0,0},
					{0,-1,1},
					{0,0,0}
			};
			diff_f_x.kernel_stencil_setup(d_f_x_kernel, 1.0/h[0]);

			double d_f_y_kernel[3][3] = {
					{0,1,0},
					{0,-1,0},
					{0,0,0},
			};
			diff_f_y.kernel_stencil_setup(d_f_y_kernel, 1.0/h[1]);


			double d_b_x_kernel[3][3] = {
					{0,0,0},
					{-1,1,0},
					{0,0,0}
			};
			diff_b_x.kernel_stencil_setup(d_b_x_kernel, 1.0/h[0]);

			double d_b_y_kernel[3][3] = {
					{0,0,0},
					{0,1,0},
					{0,-1,0},
			};
			diff_b_y.kernel_stencil_setup(d_b_y_kernel, 1.0/h[1]);


			double diff2_x_kernel[3][3] = {
					{0,0,0},
					{1.0,-2.0,1.0},
					{0,0,0}
				};
			diff2_c_x.kernel_stencil_setup(diff2_x_kernel, 1.0/(h[0]*h[0]));

			double diff2_y_kernel[3][3] = {
					{0,1.0,0},
					{0,-2.0,0},
					{0,1.0,0}
			};
			diff2_c_y.kernel_stencil_setup(diff2_y_kernel, 1.0/(h[1]*h[1]));
		}

	}

	PlaneOperators(
		PlaneDataConfig *i_planeDataConfig,		///< data config setup for spectral transformations
		const double i_domain_size[2],			///< domain size
		bool i_use_spectral_basis_diffs = true	///< use spectral differentiation (d/dx e^ix)
	)	:
		planeDataConfig(i_planeDataConfig),

		diff_c_x(i_planeDataConfig),
		diff_c_y(i_planeDataConfig),

		diff_f_x(i_planeDataConfig),
		diff_f_y(i_planeDataConfig),
		diff_b_x(i_planeDataConfig),
		diff_b_y(i_planeDataConfig),

		diff2_c_x(i_planeDataConfig),
		diff2_c_y(i_planeDataConfig),

		avg_f_x(i_planeDataConfig),
		avg_f_y(i_planeDataConfig),
		avg_b_x(i_planeDataConfig),
		avg_b_y(i_planeDataConfig),

		shift_left(i_planeDataConfig),
		shift_right(i_planeDataConfig),
		shift_up(i_planeDataConfig),
		shift_down(i_planeDataConfig)
	{
		setup(i_domain_size, i_use_spectral_basis_diffs);
	}
};



#endif /* SRC_INCLUDE_SWEET_PLANEOPERATORS_HPP_ */
