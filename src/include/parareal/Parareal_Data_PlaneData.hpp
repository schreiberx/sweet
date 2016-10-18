/*
 * PararealData_PlaneDatas.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_DATA_PLANEDATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_DATA_PLANEDATA_HPP_

#include <assert.h>
#include <parareal/Parareal_Data.hpp>



template <int N>
class Parareal_Data_PlaneData	:
		public Parareal_Data
{
public:
	PlaneData* data_arrays[N];

	Parareal_Data_PlaneData()
	{

	}


	Parareal_Data_PlaneData(
			PlaneData* i_data_arrays[N]
	)
	{
		setup(i_data_arrays);
	}


	/**
	 * Setup data
	 */
	void setup(
			PlaneData* i_data_arrays[N]
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

	const Parareal_Data&
	operator=(const Parareal_Data &i_data)
	{
		for (int i = 0; i < N; i++)
		{
			PlaneData** i_data_arrays = ((Parareal_Data_PlaneData&)i_data).data_arrays;
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

	virtual ~Parareal_Data_PlaneData()
	{
	}
};




#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_DATA_PLANEDATA_HPP_ */
