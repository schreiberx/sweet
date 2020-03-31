/*
 * SphereData.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SWEET_SPHERE_DATA_DEBUG_CONTAINER_HPP_
#define SWEET_SPHERE_DATA_DEBUG_CONTAINER_HPP_


#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>

class SphereData_DebugContainer
{
public:
	class DataContainer
	{
	public:
		std::string description;
		bool is_spectral;

		SphereData_Spectral data_spectral;
		SphereData_Physical data_physical;
	};



	static
	std::vector<DataContainer>& container_data()
	{
		static std::vector<DataContainer> foo;
		return foo;
	}


	static
	std::size_t size()
	{
		return container_data().size();
	}


	static
	void clear()
	{
		return container_data().clear();

	}


	static
	void set(
			std::size_t i,
			const SphereData_Spectral &i_data_spectral,
			const std::string &i_description
	)
	{
		if (i >= size())
			container_data().resize(i+1);

		container_data()[i].description = i_description;
		container_data()[i].is_spectral = true;
		container_data()[i].data_spectral = i_data_spectral;
	}


	static
	void set(
			std::size_t i,
			const SphereData_Physical &i_data_physical,
			const std::string &i_description
	)
	{
		if (i >= size())
			container_data().resize(i+1);

		container_data()[i].description = i_description;
		container_data()[i].is_spectral = false;
		container_data()[i].data_physical= i_data_physical;
	}


	static
	void append(
			const SphereData_Spectral &i_data_spectral,
			const std::string &i_description
	)
	{
		set(container_data().size(), i_data_spectral, i_description);
	}


	static
	void append(
			const SphereData_Physical &i_data_physical,
			const std::string &i_description
	)
	{
		set(container_data().size(), i_data_physical, i_description);
	}
};

#endif
