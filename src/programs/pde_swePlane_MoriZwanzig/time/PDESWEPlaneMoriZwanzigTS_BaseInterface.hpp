/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TS_BASE_INTERFACE_HPP_
#define SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TS_BASE_INTERFACE_HPP_

#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>

#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>
#include "ShackPDESWEPlaneTimeDiscretization.hpp"
#include "../benchmarks/ShackPDESWEMoriZwanzigPlaneBenchmarks.hpp"
#include "../ShackPDESWEPlaneMoriZwanzig.hpp"


#if SWEET_PARAREAL || SWEET_XBRAID
#include <sweet/parareal/Parareal_GenericData.hpp>
#endif

class PDESWEPlaneMoriZwanzigTS_BaseInterface
{
public:
	sweet::ErrorBase error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */

	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackExpIntegration *shackExpIntegration;
	ShackPDESWEPlaneTimeDiscretization *shackPDESWETimeDisc;
	ShackPDESWEPlaneBenchmarks *shackPDESWEBenchmark;
	ShackPDESWEPlane *shackPDESWEPlane;

	sweet::PlaneOperators *ops;

	PDESWEPlaneTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackPlaneDataOps(nullptr),
		shackIOData(nullptr),
		shackExpIntegration(nullptr),
		shackPDESWETimeDisc(nullptr),
		shackPDESWEBenchmark(nullptr),
		shackPDESWEPlane(nullptr),
		ops(nullptr)
	{

	}

	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackPlaneDataOps = io_shackDict->getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackIOData = io_shackDict->getAutoRegistration<sweet::ShackIOData>();
		shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ShackExpIntegration>();
		shackPDESWETimeDisc = io_shackDict->getAutoRegistration<ShackPDESWEPlaneTimeDiscretization>();
		shackPDESWEBenchmark = io_shackDict->getAutoRegistration<ShackPDESWEPlaneBenchmarks>();
		shackPDESWEPlane = io_shackDict->getAutoRegistration<ShackPDESWEPlane>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}


	virtual bool setup(
		sweet::PlaneOperators *io_ops
	)
	{
		ops = io_ops;
		return true;
	}


public:
	virtual void runTimestep(
			sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables

			double i_dt,		///< time step size
			double i_sim_timestamp
	) = 0;

	virtual ~PDESWEPlaneMoriZwanzigTS_BaseInterface() {}
};

#endif
