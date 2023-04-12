/*
 * PDESWEPlaneMoriZwanzig_Projection.hpp
 *
 *  Created on: 12 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_MOR_ZWANZIG_PROJECTION_HPP_
#define SRC_PROGRAMS_SWE_PLANE_MOR_ZWANZIG_PROJECTION_HPP_

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

	bool shackRegistration(sweet::ShackDictionary &io_dict)
	{
		shackPlaneDataOps = io_dict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackTimestepControl = io_dict.getAutoRegistration<sweet::ShackTimestepControl>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackIOData = io_dict.getAutoRegistration<sweet::ShackIOData>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackPDESWEPlaneMZ = io_dict.getAutoRegistration<ShackPDESWEPlane>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		return true;
	}

	void projectS(
			sweet::PlaneData_Spectral io_h_pert
			sweet::PlaneData_Spectral io_u
			sweet::PlaneData_Spectral io_v
	)
	{
	}

	void projectF(
			sweet::PlaneData_Spectral io_h_pert
			sweet::PlaneData_Spectral io_u
			sweet::PlaneData_Spectral io_v
	)
	{
	}

};
