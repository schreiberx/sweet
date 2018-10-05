/*
 * PararealData_Scalar.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_DATA_SCALAR_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_DATA_SCALAR_HPP_

#include <assert.h>
#include <parareal/Parareal_Data.hpp>



class Parareal_Data_Scalar	:
		public Parareal_Data
{
public:
	double data;

	Parareal_Data_Scalar()	:
		data(0)
	{
	}


	Parareal_Data_Scalar(
			double i_data
	)	:
		data(i_data)
	{
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


	/**
	 * Send data to rank
	 */
	void send(
			int i_mpi_rank,
			int i_mpi_comm
	)
	{
		{
//			MPI_Isend(...)
		}

		std::cout << "TODO: implement me!" << std::endl;
		assert(false);
	}

	const Parareal_Data&
	operator=(const Parareal_Data &i_data)
	{
		data = ((Parareal_Data_Scalar&)i_data).data;

		return *this;
	}

	/**
	 * Receive data from rank
	 */
	void recv(
			int i_mpi_rank,
			int i_mpi_comm
	)
	{
		{
//				MPI_Irecv(...)
		}

		std::cout << "TODO: implement me!" << std::endl;
		assert(false);
	}

	virtual ~Parareal_Data_Scalar()
	{
	}
};




#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_DATA_DATAARRAYS_HPP_ */
