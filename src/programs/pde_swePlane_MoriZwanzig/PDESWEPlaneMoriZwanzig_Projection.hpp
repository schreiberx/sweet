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
#include "PDESWEPlaneMoriZwanzig_NormalModes.hpp"
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

	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackIOData *shackIOData;

	ShackPDESWEPlaneMoriZwanzig *shackPDESWEPlaneMZ;


	double F;

	PDESWEPlaneMoriZwanzigNormalModes normal_modes;


	// min and max wavenumber modes for each wave type
	// modes_SP[0] = [min_geostr, max_geostr]   --> min_geostr <= k < max_geostr
	// modes_SP[1] = [min_west, max_west]
	// modes_SP[2] = [min_east, max_east]
	int modes_SP[3][2];
	int modes_SQ[3][2];
	int modes_SF[3][2];

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

		return true;
	}

	bool setup()
	{

		F = shackPDESWEPlaneMZ->F;

		return True;
	}


	void project_SP(
			sweet::PlaneData_Spectral io_h_pert,
			sweet::PlaneData_Spectral io_u,
			sweet::PlaneData_Spectral io_v
	)
	{

		complex eigenvalues[3];
		complex eigenvectors[3][3];


		sweet::PlaneData_Spectral h_copy = io_h_pert;
		sweet::PlaneData_Spectral u_copy = io_u_pert;
		sweet::PlaneData_Spectral v_copy = io_v_pert;

		io_h_pert.set_spectral_zero();
		io_u.set_spectral_zero();
		io_v.set_spectral_zero();

		for (int k1 = kmin; k1 < kmax; k1++)
		{
			for (int k2 = kmin; k2 < kmax; k2++)
			{
				complex U_proj[3] = {0., 0., 0.};
				normal_modes.eigendecomposition(k1, k2, eigenvalues, eigenvector);

				for (int wave_type = 0; wave_type < 3; wave_type++)
				{

					kmin = this->modes_SP[wave_type][0];
					kmax = this->modes_SP[wave_type][1];

					if (k1 < kmin || k1 > kmax || k2 < kmin || k2 > kmax)
						continue;

					complex coef_prog = 0.;
					coef_proj += eigenvector[0][wave_type] * h_copy.spectral_get(k1, k0);
					coef_proj += eigenvector[1][wave_type] * u_copy.spectral_get(k1, k0);
					coef_proj += eigenvector[2][wave_type] * v_copy.spectral_get(k1, k0);

					for (int j = 0; j < 3; j++)
						U_proj[j] += coef_proj * eigenvector[j][wave_type];

				}

				io_h_pert.spectral_set(k1, k0, U_proj[0]);
				io_u.spectral_set(k1, k0, U_proj[1]);
				io_v.spectral_set(k1, k0, U_proj[2]);
			}
		}
	}

	void project_SQ(
			sweet::PlaneData_Spectral io_h_pert,
			sweet::PlaneData_Spectral io_u,
			sweet::PlaneData_Spectral io_v
	)
	{
	}

	void project_FQ(
			sweet::PlaneData_Spectral io_h_pert,
			sweet::PlaneData_Spectral io_u,
			sweet::PlaneData_Spectral io_v
	)
	{
	}

};

#endif

