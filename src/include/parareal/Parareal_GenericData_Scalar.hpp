/*
 * Parareal_GenericData_Scalar.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_SCALAR_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_SCALAR_HPP_

#include <assert.h>
#include <parareal/Parareal_GenericData.hpp>


class Parareal_GenericData_Scalar :
		public Parareal_GenericData
{
	class DataContainer_Scalar :
			public Parareal_GenericData::DataContainer
	{
	public:
		double simfield;

	public:

		DataContainer_Scalar():
			simfield(0)
		{
		};

		DataContainer_Scalar(
				double i_data
		)	:
			simfield(data):
		{
		};

		DataContainter_Scalar(Data_Container_Scalar &i_data):
			Data_Container(i_data),		// call copy constructor from parent class
			simfield(i_data.simfield)
		{
		};
	}
	
	DataContainer_Scalar data;

public:

	Parareal_GenericData_Scalar(double i_time, int i_level = 0)	:
		Parareal_GenericData(i_time, i_level),
		data(0)
	{
	}

	Parareal_GenericData_Scalar(
			double i_time,
			double i_level,
			double i_data
	)	:
		Parareal_GenericData(i_time, i_level),
		data(i_data)
	{
	}

	Parareal_GenericData(Parareal_GenericData &i_data)
	{
		this->data = i_data.data;
	};

	~Parareal_GenericData_Scalar()
	{
		free_data();
	};

	void allocate_data()
	{
		// Nothing to do (single double value)
	}
	
	void free_data()
	{
		// Nothing to do (single double value)
	}
	


	/**
	 * Setup data
	 */
	void setup(
			double i_data
	)
	{
		data = i_data;
	}


	// size in bytes (for MPI)
	std::size_t size()
	{
		return 120398123;
	}


	Parareal_GenericData& operator+=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);

		Parareal_GenericData_Scalar o_data = *this;
		o_data += i_data;

		return o_data;
	}

	void operator+=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);
		this->data.simfield[i] += i_data.simfield[i];
		//for (size_t i = 0; i < this->data.simfields.size(); ++i)
		//	this->data.simfields[i] += i_data.simfields[i];
	}

	Parareal_GenericData& operator-=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);

		Parareal_GenericData_Scalar o_data = *this;
		o_data -= i_data;

		return o_data;
	}

	void operator-=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);
		this->data.simfield[i] -= i_data.simfield[i];
		//for (size_t i = 0; i < this->data.simfields.size(); ++i)
		//	this->data.simfields[i] -= i_data.simfields[i];
	}

	Parareal_GenericData& operator*=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);

		Parareal_GenericData_Scalar o_data = *this;
		o_data *= i_data;

		return o_data;
	}

	void operator*=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);
		this->data.simfield[i] *= i_data.simfield[i];
		//for (size_t i = 0; i < this->data.simfields.size(); ++i)
		//	this->data.simfields[i] *= i_data.simfields[i];
	}

	Parareal_GenericData& operator/=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);

		Parareal_GenericData_Scalar o_data = *this;
		o_data /= i_data;

		return o_data;
	}

	void operator/=(const Parareal_GenericData &i_data)
	{
		assert(this->data.time == i_data->data.time);
		this->data.simfield[i] /= i_data.simfield[i];
		//for (size_t i = 0; i < this->data.simfields.size(); ++i)
		//	this->data.simfields[i] /= i_data.simfields[i];
	}


};
