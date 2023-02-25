/*
 *  Created on: Dec 30, 2020
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_TRANSFORMATIONPLANS_HPP_
#define INCLUDE_SWEET_TRANSFORMATIONPLANS_HPP_

#include <string>
#include <sweet/core/ErrorBase.hpp>

namespace sweet
{

/**
 * This is a class which helps to store and load transformation plans for
 *
 *  - PlaneData_Config
 *  - SphereData_Config
 */
class TransformationPlans
{
public:
	ErrorBase error;

	enum TRANSFORMATION_PLAN_CACHE
	{
		QUICK = 1 << 1,			// quickly generate plans without caching. This is the default option.

		LOAD = 1 << 2,					// load plan if exists, otherwise generate plan
		SAVE = 1 << 3,					// save plan

		REQUIRE_LOAD = (1 << 4) | LOAD,	// force loading plan. Abort if it doesn't exist
	};


	bool getEnumFromString(const std::string &i_value, TRANSFORMATION_PLAN_CACHE& o_plan_cache)
	{
		if (i_value == "quick" || i_value == "-1")
		{
			o_plan_cache = QUICK;
			return true;
		}
		else if (i_value == "save" || i_value == "0")
		{
			o_plan_cache = TransformationPlans::SAVE;
			return true;
		}
		else if (i_value == "load" || i_value == "1")
		{
			o_plan_cache = TransformationPlans::LOAD;
			return true;
		}
		else if (i_value == "require_load" || i_value == "2")
		{
			o_plan_cache = TransformationPlans::REQUIRE_LOAD;
			return true;
		}
		else if (i_value == "load_save")
		{
			o_plan_cache = (TRANSFORMATION_PLAN_CACHE)(TransformationPlans::LOAD | TransformationPlans::SAVE);
			return true;
		}

		return error.set("Unknown option for reuse_spectral_transformation_plans '"+i_value+"'");
	}


	bool getStringFromEnum(TRANSFORMATION_PLAN_CACHE i_enum, std::string& o_plan_cache)
	{
		if (i_enum == TransformationPlans::QUICK)
		{
			o_plan_cache = "quick";
			return true;
		}
		if (i_enum == TransformationPlans::LOAD)
		{
			o_plan_cache = "load";
			return true;
		}
		if (i_enum == TransformationPlans::SAVE)
		{
			o_plan_cache = "save";
			return true;
		}
		if (i_enum == (TransformationPlans::LOAD | TransformationPlans::SAVE))
		{
			o_plan_cache = "load_save";
			return true;
		}
		if (i_enum == TransformationPlans::REQUIRE_LOAD)
		{
			o_plan_cache = "require_load";
			return true;
		}

		return error.set("Unknown enum for reuse_spectral_transformation_plans");
	}
};

}

#endif
