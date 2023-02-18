/*
 * TransformationPlans.hpp
 *
 *  Created on: Dec 30, 2020
 *      Author: Martin Schreiber
 */

#ifndef INCLUDE_SWEET_TRANSFORMATIONPLANS_HPP_
#define INCLUDE_SWEET_TRANSFORMATIONPLANS_HPP_

#include <string>


/**
 * This is a class which helps to store and load transformation plans for
 *
 *  - PlaneData_Config
 *  - SphereData_Config
 */
class TransformationPlans
{
public:
	enum TRANSFORMATION_PLAN_CACHE
	{
		QUICK = 1 << 1,			// quickly generate plans without caching. This is the default option.

		LOAD = 1 << 2,					// load plan if exists, otherwise generate plan
		SAVE = 1 << 3,					// save plan

		REQUIRE_LOAD = (1 << 4) | LOAD,	// force loading plan. Abort if it doesn't exist
	};


	static
	TRANSFORMATION_PLAN_CACHE getEnumFromString(const std::string &i_value)
	{
		if (i_value == "quick" || i_value == "-1")
			return QUICK;

		else if (i_value == "save" || i_value == "0")
			return TransformationPlans::SAVE;

		else if (i_value == "load" || i_value == "1")
			return TransformationPlans::LOAD;

		else if (i_value == "require_load" || i_value == "2")
			return TransformationPlans::REQUIRE_LOAD;

		else if (i_value == "load_save")
			return (TRANSFORMATION_PLAN_CACHE)(TransformationPlans::LOAD | TransformationPlans::SAVE);

		else
			throw std::runtime_error(std::string("Unknown option for reuse_spectral_transformation_plans '")+i_value+"'");
	}



	static
	std::string getStringFromEnum(TRANSFORMATION_PLAN_CACHE i_enum)
	{
		if (i_enum == TransformationPlans::QUICK)
			return "quick";

		if (i_enum == TransformationPlans::LOAD)
			return "load";

		if (i_enum == TransformationPlans::SAVE)
			return "save";

		if (i_enum == (TransformationPlans::LOAD | TransformationPlans::SAVE))
			return "load_save";

		if (i_enum == TransformationPlans::REQUIRE_LOAD)
			return "require_load";


		throw std::runtime_error(std::string("Unknown enum for reuse_spectral_transformation_plans"));
	}
};


#endif /* INCLUDE_SWEET_TRANSFORMATIONPLANS_HPP_ */
