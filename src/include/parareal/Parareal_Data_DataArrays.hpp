/*
 * PararealData_DataArrays.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_DATA_DATAARRAYS_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_DATA_DATAARRAYS_HPP_

#include <assert.h>
#include <parareal/Parareal_Data.hpp>



template <int N>
class PararealData_DataArrays	:
		public PararealData
{
public:
	DataArray<2>* data_arrays[N];

	PararealData_DataArrays()
	{

	}


	PararealData_DataArrays(
			DataArray<2>* i_data_arrays[N]
	)
	{
		setup(i_data_arrays);
	}


	/**
	 * Setup data
	 */
	void setup(
			DataArray<2>* i_data_arrays[N]
	)
	{
		for (int i = 0; i < N; i++)
			data_arrays[i] = i_data_arrays[i];
	}


	/**
	 * Send data to rank
	 */
	void send(
			int i_mpi_rank,
			int i_mpi_comm
	)
	{
		for (int i = 0; i < N; i++)
		{
//			MPI_Isend(...)
		}

		std::cout << "TODO: implement me!" << std::endl;
		assert(false);
	}

	const PararealData&
	operator=(const PararealData &i_data)
	{
		for (int i = 0; i < N; i++)
		{
			DataArray<2>** i_data_arrays = ((PararealData_DataArrays&)i_data).data_arrays;
			data_arrays[i] = i_data_arrays[i];
		}

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
		for (int i = 0; i < N; i++)
		{
//				MPI_Irecv(...)
		}

		std::cout << "TODO: implement me!" << std::endl;
		assert(false);
	}

	virtual ~PararealData_DataArrays()
	{
	}
};




#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_DATA_DATAARRAYS_HPP_ */
