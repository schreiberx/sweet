/*
 * PDESWEPlaneMoriZwanzig_Projection.hpp
 *
 *  Created on: 12 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_PROJECTION_HPP_
#define SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_PROJECTION_HPP_

#include <sweet/core/ErrorBase.hpp>
#include <sweet/expIntegration/ExpFunctions.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <functional>

#if SWEET_EIGEN
#include <Eigen/Eigenvalues>
#endif


#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>

#include <sweet/core/plane/PlaneSWENormalModes.hpp>
//#include "PDESWEPlaneMoriZwanzig_NormalModes.hpp"
#include "ShackPDESWEPlaneMoriZwanzig.hpp"

/**
 * SWE Plane normal mode
 */
class PDESWEPlaneMoriZwanzigProjection
{
public:

	sweet::ErrorBase error;

	typedef double T;
	typedef std::complex<T> complex;

	sweet::PlaneData_Config* planeDataConfig;

	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackIOData *shackIOData;

	ShackPDESWEPlaneMoriZwanzig *shackPDESWEPlaneMZ;

	///PDESWEPlaneMoriZwanzigNormalModes normal_modes;
	SWE_Plane_NormalModes normal_modes;


	// min and max wavenumber modes for each wave type
	// modes[j][i][0] = [min_geostr, max_geostr]   --> min_geostr <= k < max_geostr
	// modes[j][i][1] = [min_west, max_west]
	// modes[j][i][2] = [min_east, max_east]
	// i = 0 : S
	// i = 1 : F
	// j = 0: direction 0 (half the spectral size + 1)
	// j = 1: direction 1 (full spectral size)
	int modes[2][2][3][2];

	bool shackRegistration(sweet::ShackDictionary &io_dict)
	{
		shackPlaneDataOps = io_dict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackTimestepControl = io_dict.getAutoRegistration<sweet::ShackTimestepControl>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackIOData = io_dict.getAutoRegistration<sweet::ShackIOData>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackPDESWEPlaneMZ = io_dict.getAutoRegistration<ShackPDESWEPlaneMoriZwanzig>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		////normal_modes.shackRegistration(io_dict);

		return true;
	}


	bool setup(sweet::PlaneData_Config* i_planeDataConfig)
	{

		this->planeDataConfig = i_planeDataConfig;

		this->normal_modes.setup(
						shackPDESWEPlaneMZ->plane_rotating_f0,
						shackPDESWEPlaneMZ->h0,
						shackPDESWEPlaneMZ->gravitation,
						shackPlaneDataOps->plane_domain_size[0],
						shackPlaneDataOps->plane_domain_size[1],
						i_planeDataConfig
					);

		i_planeDataConfig->printInformation();
		// define S, F

		// pre-treat negative values to avoid issues with divisions by 2 below
		if (this->shackPDESWEPlaneMZ->S_geostrophic_min < 0)
			this->shackPDESWEPlaneMZ->S_geostrophic_min = -10;
		if (this->shackPDESWEPlaneMZ->S_geostrophic_max < 0)
			this->shackPDESWEPlaneMZ->S_geostrophic_max = -10;
		if (this->shackPDESWEPlaneMZ->S_gravity_west_min < 0)
			this->shackPDESWEPlaneMZ->S_gravity_west_min = -10;
		if (this->shackPDESWEPlaneMZ->S_gravity_west_max < 0)
			this->shackPDESWEPlaneMZ->S_gravity_west_max = -10;
		if (this->shackPDESWEPlaneMZ->S_gravity_east_min < 0)
			this->shackPDESWEPlaneMZ->S_gravity_east_min = -10;
		if (this->shackPDESWEPlaneMZ->S_gravity_east_max < 0)
			this->shackPDESWEPlaneMZ->S_gravity_east_max = -10;

		if (this->shackPDESWEPlaneMZ->F_geostrophic_min < 0)
			this->shackPDESWEPlaneMZ->F_geostrophic_min = -10;
		if (this->shackPDESWEPlaneMZ->F_geostrophic_max < 0)
			this->shackPDESWEPlaneMZ->F_geostrophic_max = -10;
		if (this->shackPDESWEPlaneMZ->F_gravity_west_min < 0)
			this->shackPDESWEPlaneMZ->F_gravity_west_min = -10;
		if (this->shackPDESWEPlaneMZ->F_gravity_west_max < 0)
			this->shackPDESWEPlaneMZ->F_gravity_west_max = -10;
		if (this->shackPDESWEPlaneMZ->F_gravity_east_min < 0)
			this->shackPDESWEPlaneMZ->F_gravity_east_min = -10;
		if (this->shackPDESWEPlaneMZ->F_gravity_east_max < 0)
			this->shackPDESWEPlaneMZ->F_gravity_east_max = -10;

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING

		// set physical_size = (spectral_size * 3 + 1 ) / 2  then
		// set spectral_size = physical_size

		this->modes[1][0][0][0] = (this->shackPDESWEPlaneMZ->S_geostrophic_min  * 3 + 1 ) / 2;
		this->modes[1][0][0][1] = (this->shackPDESWEPlaneMZ->S_geostrophic_max  * 3 + 1 ) / 2;
		this->modes[1][0][1][0] = (this->shackPDESWEPlaneMZ->S_gravity_west_min * 3 + 1 ) / 2;
		this->modes[1][0][1][1] = (this->shackPDESWEPlaneMZ->S_gravity_west_max * 3 + 1 ) / 2;
		this->modes[1][0][2][0] = (this->shackPDESWEPlaneMZ->S_gravity_east_min * 3 + 1 ) / 2;
		this->modes[1][0][2][1] = (this->shackPDESWEPlaneMZ->S_gravity_east_max * 3 + 1 ) / 2;

		this->modes[1][1][0][0] = (this->shackPDESWEPlaneMZ->F_geostrophic_min  * 3 + 1 ) / 2 ;
		this->modes[1][1][0][1] = (this->shackPDESWEPlaneMZ->F_geostrophic_max  * 3 + 1 ) / 2 ;
		this->modes[1][1][1][0] = (this->shackPDESWEPlaneMZ->F_gravity_west_min * 3 + 1 ) / 2 ;
		this->modes[1][1][1][1] = (this->shackPDESWEPlaneMZ->F_gravity_west_max * 3 + 1 ) / 2 ;
		this->modes[1][1][2][0] = (this->shackPDESWEPlaneMZ->F_gravity_east_min * 3 + 1 ) / 2 ;
		this->modes[1][1][2][1] = (this->shackPDESWEPlaneMZ->F_gravity_east_max * 3 + 1 ) / 2 ;
#else
		this->modes[1][0][0][0] = this->shackPDESWEPlaneMZ->S_geostrophic_min;
		this->modes[1][0][0][1] = this->shackPDESWEPlaneMZ->S_geostrophic_max;
		this->modes[1][0][1][0] = this->shackPDESWEPlaneMZ->S_gravity_west_min;
		this->modes[1][0][1][1] = this->shackPDESWEPlaneMZ->S_gravity_west_max;
		this->modes[1][0][2][0] = this->shackPDESWEPlaneMZ->S_gravity_east_min;
		this->modes[1][0][2][1] = this->shackPDESWEPlaneMZ->S_gravity_east_max;

		this->modes[1][1][0][0] = this->shackPDESWEPlaneMZ->F_geostrophic_min;
		this->modes[1][1][0][1] = this->shackPDESWEPlaneMZ->F_geostrophic_max;
		this->modes[1][1][1][0] = this->shackPDESWEPlaneMZ->F_gravity_west_min;
		this->modes[1][1][1][1] = this->shackPDESWEPlaneMZ->F_gravity_west_max;
		this->modes[1][1][2][0] = this->shackPDESWEPlaneMZ->F_gravity_east_min;
		this->modes[1][1][2][1] = this->shackPDESWEPlaneMZ->F_gravity_east_max;

#endif


		// treat negative values
		for (int i = 0; i < 2; i++)
			for (int wave_type = 0; wave_type < 3; wave_type++)
				for (int j = 0; j < 2; j++)
					if (this->modes[1][i][wave_type][j] < 0)
						this->modes[1][i][wave_type][j] = 0;

		// modes in direction 1 (half the spectral size + 1)
		for (int i = 0; i < 2; i++)
			for (int wave_type = 0; wave_type < 3; wave_type++)
				for (int j = 0; j < 2; j++)
				{
					if (this->modes[1][i][wave_type][j] == 0)
						this->modes[0][i][wave_type][j] = 0;
					else
						this->modes[0][i][wave_type][j] = this->modes[1][i][wave_type][j] / 2 + 1;
				}

		std::vector<std::string> strs = {"S", "F"};
		std::vector<std::string> strw = {"geostrophic", "gravity west", "gravity east"};

		std::cout << "DEFINITION OF S AND F" << std::endl;
		for (int i = 0; i < 2; i++)
			for (int wave_type = 0; wave_type < 3; wave_type++)
				std::cout << strs[i] << ", " << strw[wave_type] << ": " << modes[0][i][wave_type][0] << " <= k1 < " << modes[0][i][wave_type][1] << "; " << modes[1][i][wave_type][0] << " <= k2 < " << modes[1][i][wave_type][1] << std::endl;


		// Check if S and F overlap
		for (int i = 0; i < 2; i++)
		{
			for (int j = i + 1; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					for (int wave_type = 0; wave_type < 3; wave_type++)
					{
						int min1 = this->modes[k][i][wave_type][0];
						int max1 = this->modes[k][i][wave_type][1];
						int min2 = this->modes[k][j][wave_type][0];
						int max2 = this->modes[k][j][wave_type][1];

						if (
							(min2 < max1 && max2 > max1) ||
							(min2 < min1 && max2 > min1) ||
							(min1 < max2 && max1 > max2) ||
							(min1 < max2 && max1 > max2)
						)
							SWEETError("Regions " + strs[i] + " and " + strs[j] + " overlap in " + strw[wave_type] + " waves in direction " + std::to_string(k) + " !");
					}
				}
			}
		}

		// check if SP, SQ, FQ cover the entire spectral space
		for (int wave_type = 0; wave_type < 3; wave_type++)
		{
			for (int k = 0; k < 2; k++)
			{
				// check if there is a min = zero
				bool zero_ok = false;
				for (int i = 0; i < 2; i++)
					if (this->modes[k][i][wave_type][0] == 0 && this->modes[k][i][wave_type][1] > 0)
					{
						zero_ok = true;
						break;
					}
				if (!zero_ok)
					SWEETError(strw[wave_type] + " waves are not fully covered by S and F in direction " + std::to_string(k) + " !");

				// sum modes to check if they are all taken into account
				int sum = 0;
				for (int i = 0; i < 2; i++)
				{
					int min = this->modes[k][i][wave_type][0];
					int max = this->modes[k][i][wave_type][1] - 1;
					sum += (min + max) * (max - min + 1) / 2;
				}
				if (sum != (i_planeDataConfig->spectral_data_size[k] - 1) * i_planeDataConfig->spectral_data_size[k] / 2)
					SWEETError(strw[wave_type] + " waves are not fully covered by S and F in direction " + std::to_string(k) + " !");
			}
		}



		return true;
	}



	void project(
			sweet::PlaneData_Spectral &io_h_pert,
			sweet::PlaneData_Spectral &io_u,
			sweet::PlaneData_Spectral &io_v,
			std::string projection_type
	)
	{

		int idx;
		if (projection_type == "S")
			idx = 0;
		else if (projection_type == "F")
			idx = 1;
		else
			SWEETError("Wrong projection type: " + projection_type );

		complex eigenvectors[3][3];

		sweet::PlaneData_Spectral h_copy = io_h_pert;
		sweet::PlaneData_Spectral u_copy = io_u;
		sweet::PlaneData_Spectral v_copy = io_v;

		io_h_pert.spectral_set_zero();
		io_u.spectral_set_zero();
		io_v.spectral_set_zero();

		// get min and max K from all wave types (GEO, GW, GE)
		int Kmin1 = std::min(std::min(this->modes[0][idx][0][0], this->modes[0][idx][1][0]), this->modes[0][idx][2][0]);
		int Kmax1 = std::max(std::max(this->modes[0][idx][0][1], this->modes[0][idx][1][1]), this->modes[0][idx][2][1]);
		int Kmin2 = std::min(std::min(this->modes[1][idx][0][0], this->modes[1][idx][1][0]), this->modes[1][idx][2][0]);
		int Kmax2 = std::max(std::max(this->modes[1][idx][0][1], this->modes[1][idx][1][1]), this->modes[1][idx][2][1]);

		for (int k1 = Kmin1; k1 < Kmax1; k1++)
		{

			////// Only half of the modes are defined in k1 direction!!
			////if (k1 >= this->planeDataConfig->spectral_data_size[0])
			////	continue;

			for (int k2 = Kmin2; k2 < Kmax2; k2++)
			{

				complex U_proj[3] = {0., 0., 0.};

				normal_modes.eigendecomposition(k1, k2, eigenvectors);

				///////check if eigenvectors are orthonormal
				/////for (int i = 0; i < 3; i++)
				/////	for (int j = 0; j < 3; j++)
				/////	{
				/////		complex sum = 0.;
				/////		for (int k = 0; k < 3; k++)
				/////			sum += eigenvectors[k][i] * std::conj(eigenvectors[k][j]);
				/////		std::cout << "AAAAAAAAA " << k1 << " " << k2 << " " << i << " " << j << " " << sum << std::endl;
				/////	}
				/////std::cout << std::endl;

				for (int wave_type = 0; wave_type < 3 ; wave_type++)
				{

					int kmin1 = this->modes[0][idx][wave_type][0];
					int kmax1 = this->modes[0][idx][wave_type][1];
					int kmin2 = this->modes[1][idx][wave_type][0];
					int kmax2 = this->modes[1][idx][wave_type][1];

					if (k1 < kmin1 || k1 >= kmax1 || k2 < kmin2 || k2 >= kmax2)
						continue;

					complex coef_proj = 0.;
					coef_proj += std::conj(eigenvectors[0][wave_type]) * h_copy.spectral_get(k2, k1);
					coef_proj += std::conj(eigenvectors[1][wave_type]) * u_copy.spectral_get(k2, k1);
					coef_proj += std::conj(eigenvectors[2][wave_type]) * v_copy.spectral_get(k2, k1);

					for (int j = 0; j < 3; j++)
						U_proj[j] += coef_proj * eigenvectors[j][wave_type];

				}

				/////std::cout << "AAAAAA " << projection_type << " " << k2 << " " << k1 << " " << U_proj[0] << std::endl;
				io_h_pert.spectral_set(k2, k1, U_proj[0]);
				io_u.spectral_set(k2, k1, U_proj[1]);
				io_v.spectral_set(k2, k1, U_proj[2]);
			}
		}

		io_h_pert.spectral_zeroAliasingModes();
		io_u.spectral_zeroAliasingModes();
		io_v.spectral_zeroAliasingModes();

	}

	void project_S(
			sweet::PlaneData_Spectral &io_h_pert,
			sweet::PlaneData_Spectral &io_u,
			sweet::PlaneData_Spectral &io_v
	)
	{
		this->project(io_h_pert, io_u, io_v, "S");
	}

	void project_F(
			sweet::PlaneData_Spectral &io_h_pert,
			sweet::PlaneData_Spectral &io_u,
			sweet::PlaneData_Spectral &io_v
	)
	{
		this->project(io_h_pert, io_u, io_v, "F");
	}

};

#endif

