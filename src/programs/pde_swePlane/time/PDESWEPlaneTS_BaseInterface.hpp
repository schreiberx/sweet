/*
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TS_BASE_INTERFACE_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TS_BASE_INTERFACE_HPP_

#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>

#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>
#include "ShackPDESWEPlaneTimeDiscretization.hpp"
#include "../benchmarks/ShackPDESWEPlaneBenchmarks.hpp"
#include "../ShackPDESWEPlane.hpp"


#if SWEET_PARAREAL || SWEET_XBRAID
#include <sweet/parareal/Parareal_GenericData.hpp>
#endif

class PDESWEPlaneTS_BaseInterface
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
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*io_shackDict);

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
			sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,		///< time step size
			double i_sim_timestamp
	) = 0;

#if (SWEET_PARAREAL && SWEET_PARAREAL_PLANE) || (SWEET_XBRAID && SWEET_XBRAID_PLANE)
	void runTimestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		sweet::PlaneData_Spectral h_pert = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]);
		sweet::PlaneData_Spectral u = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]);
		sweet::PlaneData_Spectral v = *(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2]);

		runTimestep(h_pert, u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0]) = h_pert;
		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1]) = u;
		*(io_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2]) = v;

	}

	// for parareal SL
	virtual void set_previous_solution(
				sweet::PlaneData_Spectral &i_h_prev,
				sweet::PlaneData_Spectral &i_u_prev,
				sweet::PlaneData_Spectral &i_v_prev
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			Parareal_GenericData* i_data
	)
	{
		sweet::PlaneData_Spectral h_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[0];
		sweet::PlaneData_Spectral u_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[1];
		sweet::PlaneData_Spectral v_prev = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[2];

		set_previous_solution(h_prev, u_prev, v_prev);
	};
#endif

	virtual ~PDESWEPlaneTS_BaseInterface() {}
};

#endif
