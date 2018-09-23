/*
 * SphereDiagnostics.hpp
 *
 *  Created on: 25 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_SPHEREDIAGNOSTICS_HPP_
#define SRC_INCLUDE_SWEET_PLANE_SPHEREDIAGNOSTICS_HPP_

#include <sweet/sphere/SphereDataConfig.hpp>


class SphereDiagnostics
{
	SphereDataConfig *sphereDataConfig;
	SphereData modeIntegralValues;

	SphereDataPhysical fg;

	/*
	 * Gaussian quadrature weights
	 */
	std::vector<double> gauss_weights;


public:
	SphereDiagnostics(
			SphereDataConfig *i_sphereDataConfig,
			SimulationVariables &i_simVars,
			int i_verbose = 1
	)	:
		sphereDataConfig(i_sphereDataConfig),
		modeIntegralValues(sphereDataConfig),
		fg(sphereDataConfig)
	{
		gauss_weights.resize(sphereDataConfig->physical_num_lat);

		SphereData modeSelector(sphereDataConfig);
		modeSelector.spectral_set_zero();

		if (i_verbose > 0)
			std::cout << "Setting up SphereDiagnostics" << std::endl;

		int n = shtns_gauss_wts(sphereDataConfig->shtns, gauss_weights.data());

		if (n*2 != sphereDataConfig->physical_num_lat)
		{
			std::cerr << "Returned " << n << " number of Gaussian quadrature point (halved)" << std::endl;
			FatalError("Wrong number of Gaussian quadrature points given!");
		}

		for (int i = 0; i < sphereDataConfig->physical_num_lat/2; i++)
			gauss_weights[sphereDataConfig->physical_num_lat-i-1] = gauss_weights[i];


#if 0
		for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
		{
			if (i_verbose > 0)
				std::cout << "Mode m = " << m << " / " << sphereDataConfig->spectral_modes_m_max << std::endl;

			std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
			{
				double integral = 0;

#if !SWEET_DEBUG
				if (m != 0)
#endif
				{
					modeSelector.spectral_set_zero();

					// activate mode
					modeSelector.spectral_space_data[idx] = 1.0;

					// conversion to physical space stores Gaussian quadrature weights
					modeSelector.request_data_physical();


					integral = modeSelector.physical_reduce_sum_quad();

					if (m != 0 && integral > 10e-12)
					{
						std::cout << n << ", " << m << ": Integral value expected to be close to zero, but is " << integral << std::endl;
						FatalError("Integral value not close to zero for m != 0");
					}
				}

				modeIntegralValues.spectral_space_data[idx] = integral;
				idx++;
			}
		}
#endif

		fg.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu*2.0*i_simVars.sim.coriolis_omega;
			}
		);

		/*
		 * Test Gaussian quadrature
		 *
		 * Accurate for order (2n-1)
		 *
		 * test with y(x) = x^(2n-1)
		 */
		{
			int n = sphereDataConfig->physical_num_lat;
			double sum = 0;
			for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
			{
				double x = sphereDataConfig->lat_gaussian[jlat];
				double value = std::pow(x, 2.0*n-1.0) + std::pow(x, 2.0*n-2.0);

				sum += value*gauss_weights[jlat];
			}

			double accurate = 0;
			accurate += 1.0/(2.0*n)*(std::pow(1.0, 2.0*n)-std::pow(-1.0, 2.0*n));
			accurate += 1.0/(2.0*n-1.0)*(std::pow(1.0, 2.0*n-1.0)-std::pow(-1.0, 2.0*n-1.0));

			if (std::abs(accurate-sum) > 1e-10)
				FatalError("Error in quadrature");
		}
	}



#if 0

public:
	double compute_sphere_integral(
			const SphereDataPhysical &i_data
	)	const
	{
		double sum = 0;

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
#error "TODO"
#else

		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double value = i_data.physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon];

				value *= gauss_weights[jlat]*sphereDataConfig->lat_cogaussian[jlat];

				sum += value;
			}
		}
#endif

		sum /= (double)sphereDataConfig->physical_num_lon;

		sum *= (4.0/M_PI);

		return sum;
	}

#endif


public:
	/*
	 * The integral similar to the zylinder is used because of the lat-related scaling factor.
	 */
	double compute_zylinder_integral(
			const SphereData &i_data
	)	const
	{
		i_data.request_data_physical();

		double sum = 0;

#if SPHERE_DATA_GRID_LAYOUT	== SPHERE_DATA_LAT_CONTINUOUS
#error "TODO"
#else
		for (int jlat = 0; jlat < sphereDataConfig->physical_num_lat; jlat++)
		{
			for (int ilon = 0; ilon < sphereDataConfig->physical_num_lon; ilon++)
			{
				double value = i_data.physical_space_data[jlat*sphereDataConfig->physical_num_lon + ilon];

				sum += value*gauss_weights[jlat];
			}
		}
#endif
		sum /= (double)sphereDataConfig->physical_num_lon;

//		sum *= 0.5;
//		sum *= M_PI*4.0;
		sum *= 2.0*M_PI;

		return sum;
	}



public:
	void update_h_u_v_2_mass_energy_enstrophy_4_zylinder(
			const SphereOperators &op,
			const SphereData &i_prog_h,
			const SphereData &i_prog_u,
			const SphereData &i_prog_v,

			SimulationVariables &io_simVars
	)
	{
		SphereDataPhysical h = i_prog_h.getSphereDataPhysical();
		SphereDataPhysical u = i_prog_u.getSphereDataPhysical();
		SphereDataPhysical v = i_prog_v.getSphereDataPhysical();

		// Convert to non-robert formulation
		// no problem, since we do everything in physical space
		u = u.robert_convertToNonRobert();
		v = v.robert_convertToNonRobert();

		double normalization = 4.0*M_PI*(io_simVars.sim.earth_radius*io_simVars.sim.earth_radius);

		// mass
		io_simVars.diag.total_mass = compute_zylinder_integral(h) * normalization;

		// energy
		SphereDataPhysical pot_energy = h*(io_simVars.sim.gravitation*normalization);
		SphereDataPhysical kin_energy = h*(u*u+v*v)*(0.5*normalization);

		io_simVars.diag.potential_energy = compute_zylinder_integral(pot_energy);
		io_simVars.diag.kinetic_energy = compute_zylinder_integral(kin_energy);

		io_simVars.diag.total_energy = io_simVars.diag.kinetic_energy + io_simVars.diag.potential_energy;

		// total vorticity
		SphereDataPhysical eta(i_prog_h.sphereDataConfig);
//		if (io_simVars.misc.sphere_use_robert_functions)
//			eta = op.robert_uv_to_vort(u, v).getSphereDataPhysical();
//		else
			eta = op.uv_to_vort(u, v).getSphereDataPhysical();

		eta += fg;

		// enstrophy
		io_simVars.diag.total_potential_enstrophy = 0.5*compute_zylinder_integral(eta*eta) * normalization;
	}




public:
	void update_phi_vort_div_2_mass_energy_enstrophy(
			const SphereOperators &op,
			const SphereData &i_prog_phi,
			const SphereData &i_prog_vort,
			const SphereData &i_prog_div,

			SimulationVariables &io_simVars
	)
	{
		SphereDataPhysical h(sphereDataConfig);
		SphereDataPhysical u(sphereDataConfig);
		SphereDataPhysical v(sphereDataConfig);

		h = i_prog_phi.getSphereDataPhysical()*(1.0/io_simVars.sim.gravitation);
		if (io_simVars.misc.sphere_use_robert_functions)
			op.robert_vortdiv_to_uv(i_prog_vort, i_prog_div, u, v);
		else
			op.vortdiv_to_uv(i_prog_vort, i_prog_div, u, v);

		double normalization = (io_simVars.sim.earth_radius*io_simVars.sim.earth_radius);

		// mass
		io_simVars.diag.total_mass = compute_zylinder_integral(h) * normalization;

		// energy
		//SphereDataPhysical pot_energy = h*(io_simVars.sim.gravitation*normalization);
		SphereDataPhysical pot_energy = h*h*0.5*normalization;
		SphereDataPhysical kin_energy = h*(u*u+v*v)*(0.5*normalization);

		io_simVars.diag.potential_energy = compute_zylinder_integral(pot_energy);
		io_simVars.diag.kinetic_energy = compute_zylinder_integral(kin_energy);

		/*
		 * We follow the Williamson et al. equation (137) here
		 */
//		double dummy_energy = compute_zylinder_integral(h*h*(0.5*normalization));
//		io_simVars.diag.total_energy = io_simVars.diag.kinetic_energy + dummy_energy;//io_simVars.diag.potential_energy;

		/*
		 * We follow pot/kin energy here
		 */
		io_simVars.diag.total_energy = io_simVars.diag.kinetic_energy + io_simVars.diag.potential_energy;

		// total vorticity
		// TODO: maybe replace this with the i_vort parameter
		SphereDataPhysical eta(h.sphereDataConfig);
		if (io_simVars.misc.sphere_use_robert_functions)
			eta = op.robert_uv_to_vort(u, v).getSphereDataPhysical();
		else
			eta = op.uv_to_vort(u, v).getSphereDataPhysical();

		eta += fg;

		// enstrophy (Williamson paper, equation 138)
		io_simVars.diag.total_potential_enstrophy = 0.5*compute_zylinder_integral(eta*eta/h) * normalization;
	}


};



#endif /* SRC_INCLUDE_SWEET_PLANE_SPHEREDIAGNOSTICS_HPP_ */
