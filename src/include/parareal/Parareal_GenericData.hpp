/*
 * Parareal_GenericData.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_HPP_

/**
 * Generic Parareal Data class allowing a common parareal interface for any kind of data (plane, sphere);
 * This class may be inherited and specialized to each data type
 */

class Parareal_GenericData
{
	/*
	 * For Parareal data
	 * start of coarse time step
	 * end of coarse time step
	 * corrections
	 */

	class DataContainer
	{
	public:
		int level;
		double time;

	public:

		Data_Container(double i_time, int i_level = 0)
			level(i_level),
			time(i_time)
		{
		};

		Data_Containter(Data_Container &i_data) :
			level(i_data.level),
			time(i_data.time)
		{
		};
	}


public:

	Parareal_GenericData(double i_time, int i_level = 0)
	{
		this->data.time = i_time;
		this->data.level = i_level;
	};

	virtual Parareal_GenericData(Parareal_GenericData &i_data) = 0;
	virtual Parareal_GenericData(Parareal_GenericData &&i_data) = 0;

	virtual void allocate_data()=0;
	
	virtual void free_data()=0;

#if SWEET_MPI
	virtual std::size_t size() = 0;
	virtual void serialize(void *data) = 0;
	virtual void deserialize(void *data) = 0;
#endif


	virtual Parareal_GenericData& operator+(const Parareal_GenericData &i_data) = 0;
	virtual void operator+=(const Parareal_GenericData &i_data) = 0;

	virtual Parareal_GenericData& operator-(const Parareal_GenericData &i_data) = 0;
	virtual void operator-=(const Parareal_GenericData &i_data) = 0;

	virtual Parareal_GenericData& operator*(const Parareal_GenericData &i_data) = 0;
	virtual void operator*=(const Parareal_GenericData &i_data) = 0;

	virtual Parareal_GenericData& operator/(const Parareal_GenericData &i_data) = 0;
	virtual void operator/=(const Parareal_GenericData &i_data) = 0;
};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_HPP_ */
