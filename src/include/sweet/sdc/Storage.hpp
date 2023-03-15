/*
 * Storage.hpp
 *
 *  Created on: 7 March 2023
 *      Author: Thibaut LUNET <thibaut.lunet@tuhh.de>
 */
#ifndef SWEET_SDC_STORAGE_HPP_
#define SWEET_SDC_STORAGE_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>

/*
 * Class to store solution data at one node
 */
class SWE_VariableVector
{
public:
	sweet::SphereData_Spectral phi;
	sweet::SphereData_Spectral vrt;
	sweet::SphereData_Spectral div;

public:

	// Copy constructor
	SWE_VariableVector(const SWE_VariableVector &i_value)	:
		phi(i_value.phi.sphereDataConfig),
		vrt(i_value.phi.sphereDataConfig),
		div(i_value.phi.sphereDataConfig)
	{
		phi = i_value.phi;
		vrt = i_value.vrt;
		div = i_value.div;
	}

	// Default constructor
	SWE_VariableVector()
	{
	}

	// Fill values for phi, vort and div
	SWE_VariableVector& operator=(const SWE_VariableVector& u) {
		phi = u.phi;
		vrt = u.vrt;
		div = u.div;
		return *this;
	}

	bool setup(const sweet::SphereData_Config *sphere_data_config)
	{
		phi.setup(sphere_data_config);
		vrt.setup(sphere_data_config);
		div.setup(sphere_data_config);

		return true;
	}

	void swap(SWE_VariableVector &io_value)
	{
		phi.swap(io_value.phi);
		vrt.swap(io_value.vrt);
		div.swap(io_value.div);
	}
};


// Class to store all the solution data to each nodes and two iterations
class SDC_NodeStorage {
	std::vector<SWE_VariableVector> data;


public:
	SDC_NodeStorage()
	{
	}

public:
	void setup(
			const sweet::SphereData_Config* sphereDataConfig,
			size_t num_nodes
	){
		data.resize(num_nodes);

		for (size_t i = 0; i < num_nodes; i++)
		{
			data[i].setup(sphereDataConfig);
		}
	}

	SWE_VariableVector& operator[](int i)
	{
		return data[i];
	}

	void swap(SDC_NodeStorage &i_value)
	{
		for (size_t i = 0; i < data.size(); i++)
			data[i].swap(i_value.data[i]);
	}
};

#endif