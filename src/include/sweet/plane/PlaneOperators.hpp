/*
 * PlaneOperators.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_PLANEOPERATORS_HPP_
#define SRC_INCLUDE_SWEET_PLANEOPERATORS_HPP_


#if SWEET_USE_PLANE_SPECTRAL_SPACE
//	#include <sweet/plane/PlaneDataComplex.hpp>
#endif

#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneDataConfig.hpp>


class PlaneOperators
{
public:
	PlaneDataConfig *planeDataConfig;

public:
	// differential operators (central / forward / backward)
	PlaneData_Spectral diff_c_x, diff_c_y;
	PlaneData_Spectral diff_f_x, diff_f_y;
	PlaneData_Spectral diff_b_x, diff_b_y;

	PlaneData_Spectral diff2_c_x, diff2_c_y;

	PlaneData_Physical avg_f_x, avg_f_y;
	PlaneData_Physical avg_b_x, avg_b_y;

	PlaneData_Physical shift_left;
	PlaneData_Physical shift_right;
	PlaneData_Physical shift_up;
	PlaneData_Physical shift_down;


	/**
	 * D2, e.g. for viscosity
	 */
	PlaneData_Spectral diff2(
			const PlaneData_Spectral &i_dataArray
	)
	{
		return diff2_c_x(i_dataArray) + diff2_c_y(i_dataArray);
	}



	/**
	 *        __2
	 * apply  \/  operator (aka Laplace)
	 */
	inline PlaneData_Spectral laplace(
			const PlaneData_Spectral &i_a
	)
	{
		return diff2_c_x(i_a)+diff2_c_y(i_a);
	}



	/**
	 * Vorticity
	 *
	 * vort(a,b) = db/dx - da/dy
	 */
	PlaneData_Spectral vort(
			const PlaneData_Spectral &a,
			const PlaneData_Spectral &b
	)
	{
		return diff_c_x(b) - diff_c_y(a);
	}



	/**
	 * Divergence
	 *
	 * div(a,b) = da/dx + db/dy
	 */
	PlaneData_Spectral div(
			const PlaneData_Spectral &a,
			const PlaneData_Spectral &b
	)
	{
		return diff_c_x(a) + diff_c_y(b);
	}

	/**
	 * kinetic energy
	 *
	 * ke(a,b) = 0.5*(a^2+b^2)
	 */
	PlaneData_Spectral ke(
			const PlaneData_Spectral &a,
			const PlaneData_Spectral &b
	)
	{
		return 0.5*(a*a+b*b);
	}



	/*
	 * Compute Arakawa Jacobian
	 * See A. Arakawa, V. R. Lamb, "A potential enstrophy and energy conserving scheme for the shallow water equations"
	 *
	 * J(a,b) = da/dx db/dy - da/dy db/dx
	 */
	PlaneData_Spectral J(
			const PlaneData_Spectral &a,
			const PlaneData_Spectral &b
	)
	{
		return diff_c_x(a)*diff_c_y(b) - diff_c_y(a)*diff_c_x(b);
	}


	/*
	 * Compute time derivative of Arakawa Jacobian
	 *
	 * J(a,b)_t = (da/dx db/dy - da/dy db/dx)_t
	 */
	PlaneData_Spectral J_t(
			const PlaneData_Spectral &a,
			const PlaneData_Spectral &b,
			const PlaneData_Spectral &a_t,
			const PlaneData_Spectral &b_t
	)
	{
		return	  diff_c_x(a_t)*diff_c_y(b)
				+ diff_c_x(a)*diff_c_y(b_t)
				- diff_c_y(a_t)*diff_c_x(b)
				- diff_c_y(a)*diff_c_x(b_t);
	}



	/**
	 *        __
	 * apply  \/ .  operator
	 */
	inline PlaneData_Spectral diff_dot(
			const PlaneData_Spectral &i_a
	)
	{
		return diff_c_x(i_a)+diff_c_y(i_a);
	}



	/**
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	inline PlaneData_Spectral diffN_x(
			const PlaneData_Spectral &io_u,
			int i_order
	)
	{
		if (i_order == 0)
			return io_u;

		PlaneData_Spectral tu = io_u;

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
	inline PlaneData_Spectral diffN_y(
			const PlaneData_Spectral &io_v,
			int i_order
	)
	{
		if (i_order == 0)
			return io_v;

		PlaneData_Spectral tv = io_v;

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
	inline PlaneData_Spectral diffusion_coefficient(
			int i_order
	)
	{
		//Check if even
		assert( i_order % 2 == 0);
		assert( i_order > 0);
		PlaneData_Spectral out = diff2_c_x+diff2_c_y;

#if 0
		for (int i = 1; i < i_order/2; i++)
			out = pow(-1, i)*(diff2_c_x(out)+diff2_c_y(out));
#else
		/*
		 * Always use negative sign for hyperdiffusion to allow using always positive viscosity
		 */
		for (int i = 1; i < i_order/2; i++)
			out = -(diff2_c_x(out)+diff2_c_y(out));
#endif

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
	inline PlaneData_Spectral implicit_diffusion(
			const PlaneData_Spectral &i_data,
			double i_coef,
			int i_order
	)
	{
		PlaneData_Spectral out=i_data;

		// Get diffusion coefficients (these are the -mu*dt*D^q, where q is the order
		PlaneData_Spectral diff = -i_coef*diffusion_coefficient(i_order);

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
			diff_c_x.spectral_set_zero();

			for (int j = planeDataConfig->spectral_data_iteration_ranges[0][1][0]; j < (int)planeDataConfig->spectral_data_iteration_ranges[0][1][1]; j++)
			{
				for (int i = planeDataConfig->spectral_data_iteration_ranges[0][0][0]; i < (int)planeDataConfig->spectral_data_iteration_ranges[0][0][1]; i++)
				{
					std::complex<double> data(0.0, ((double)i*2.0*M_PI/(double)i_domain_size[0]));
					diff_c_x.spectral_set(j, i, data);
				}
			}

			for (int j = planeDataConfig->spectral_data_iteration_ranges[1][1][0]; j < (int)planeDataConfig->spectral_data_iteration_ranges[1][1][1]; j++)
			{
				for (int i = planeDataConfig->spectral_data_iteration_ranges[1][0][0]; i < (int)planeDataConfig->spectral_data_iteration_ranges[1][0][1]; i++)
				{
					std::complex<double> data(0.0, ((double)i*2.0*M_PI/(double)i_domain_size[0]));
					diff_c_x.spectral_set(j, i, data);
				}
			}

			/*
			 * DIFF operator in y axis
			 */
			diff_c_y.spectral_set_zero();

			for (int j = planeDataConfig->spectral_data_iteration_ranges[0][1][0]; j < (int)planeDataConfig->spectral_data_iteration_ranges[0][1][1]; j++)
			{
				for (int i = planeDataConfig->spectral_data_iteration_ranges[0][0][0]; i < (int)planeDataConfig->spectral_data_iteration_ranges[0][0][1]; i++)
				{
					std::complex<double> data(0, (double)((double)j*2.0*M_PI/(double)i_domain_size[1]));
					diff_c_y.spectral_set(j, i, data);
				}
			}


			for (int j = planeDataConfig->spectral_data_iteration_ranges[1][1][0]; j < (int)planeDataConfig->spectral_data_iteration_ranges[1][1][1]; j++)
			{
				for (int i = planeDataConfig->spectral_data_iteration_ranges[1][0][0]; i < (int)planeDataConfig->spectral_data_iteration_ranges[1][0][1]; i++)
				{
					std::complex<double> data(0, -(double)((double)(planeDataConfig->spectral_data_size[1]-j)*2.0*M_PI/(double)i_domain_size[1]));
					diff_c_y.spectral_set(j, i, data);
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
			diff_f_x.toPhys().kernel_stencil_setup(d_f_x_kernel, 1.0/h[0]);

			double d_f_y_kernel[3][3] = {
					{0,1,0},
					{0,-1,0},
					{0,0,0},
			};
			diff_f_y.toPhys().kernel_stencil_setup(d_f_y_kernel, 1.0/h[1]);


			double d_b_x_kernel[3][3] = {
					{0,0,0},
					{-1,1,0},
					{0,0,0}
			};
			diff_b_x.toPhys().kernel_stencil_setup(d_b_x_kernel, 1.0/h[0]);

			double d_b_y_kernel[3][3] = {
					{0,0,0},
					{0,1,0},
					{0,-1,0},
			};
			diff_b_y.toPhys().kernel_stencil_setup(d_b_y_kernel, 1.0/h[1]);


			/*
			 * 2nd order differential operators
			 */
			diff2_c_x.spectral_set_zero();
			for (int r = 0; r < 2; r++)
			{
				for (int j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < (int)planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
				{
					for (int i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < (int)planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
					{
						std::complex<double> data = diff_c_x.spectral_get(j, i);
						diff2_c_x.spectral_set(j, i, data*data);
					}
				}
			}

			diff2_c_y.spectral_set_zero();
			for (int r = 0; r < 2; r++)
			{
				for (int j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < (int)planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
				{
					for (int i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < (int)planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
					{
						std::complex<double> data = diff_c_y.spectral_get(j, i);
						diff2_c_y.spectral_set(j, i, data*data);
					}
				}
			}

			diff_c_x.spectral_debugCheckForZeroAliasingModes();
			diff_c_y.spectral_debugCheckForZeroAliasingModes();
			diff2_c_x.spectral_debugCheckForZeroAliasingModes();
			diff2_c_y.spectral_debugCheckForZeroAliasingModes();
#endif
		}
		else
		{
			PlaneData_Physical tmp(planeDataConfig);

			double diff1_x_kernel[3][3] = {
					{0,0,0},
					{-1.0,0,1.0},
					{0,0,0}
			};
			tmp.kernel_stencil_setup(diff1_x_kernel, 1.0/(2.0*h[0]));
			diff_c_x.loadPlaneDataPhysical(tmp);

			double diff1_y_kernel[3][3] = {
					{0,1.0,0},	// higher y coordinate
					{0,0,0},
					{0,-1.0,0},	// lower y coordinate
			};
			tmp.kernel_stencil_setup(diff1_y_kernel, 1.0/(2.0*h[1]));
			diff_c_y.loadPlaneDataPhysical(tmp);

			double d_f_x_kernel[3][3] = {
					{0,0,0},
					{0,-1,1},
					{0,0,0}
			};
			tmp.kernel_stencil_setup(d_f_x_kernel, 1.0/h[0]);
			diff_f_x.loadPlaneDataPhysical(tmp);

			double d_f_y_kernel[3][3] = {
					{0,1,0},
					{0,-1,0},
					{0,0,0},
			};
			tmp.kernel_stencil_setup(d_f_y_kernel, 1.0/h[1]);
			diff_f_y.loadPlaneDataPhysical(tmp);


			double d_b_x_kernel[3][3] = {
					{0,0,0},
					{-1,1,0},
					{0,0,0}
			};
			tmp.kernel_stencil_setup(d_b_x_kernel, 1.0/h[0]);
			diff_b_x.loadPlaneDataPhysical(tmp);

			double d_b_y_kernel[3][3] = {
					{0,0,0},
					{0,1,0},
					{0,-1,0},
			};
			tmp.kernel_stencil_setup(d_b_y_kernel, 1.0/h[1]);
			diff_b_y.loadPlaneDataPhysical(tmp);


			double diff2_x_kernel[3][3] = {
					{0,0,0},
					{1.0,-2.0,1.0},
					{0,0,0}
				};
			tmp.kernel_stencil_setup(diff2_x_kernel, 1.0/(h[0]*h[0]));
			diff2_c_x.loadPlaneDataPhysical(tmp);

			double diff2_y_kernel[3][3] = {
					{0,1.0,0},
					{0,-2.0,0},
					{0,1.0,0}
			};
			tmp.kernel_stencil_setup(diff2_y_kernel, 1.0/(h[1]*h[1]));
			diff2_c_y.loadPlaneDataPhysical(tmp);
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
