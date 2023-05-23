#ifndef SRC_PROGRAMS_SIMDATA_TIMESTEPPER_BASE_HELPER_HPP_
#define SRC_PROGRAMS_SIMDATA_TIMESTEPPER_BASE_HELPER_HPP_

#include "DESolver_DataContainer_Base.hpp"

namespace sweet
{

template <typename TData, typename TConfig>
class DESolver_TimeTreeNode_CastHelper
{
public:
	/*
	 * Casting functions which are used for convenience
	 */
	static inline
	TData& cast(sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<TData&>(i_U);
	}

	static inline
	const TData& cast(const sweet::DESolver_DataContainer_Base &i_U)
	{
		return static_cast<const TData&>(i_U);
	}

	const TConfig& cast(
			const sweet::DESolver_Config_Base& i_deTermConfig
	)
	{
		return static_cast<const TConfig&>(i_deTermConfig);
	}
};

}

#endif
