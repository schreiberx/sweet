/*
 * PararealData.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_DATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_DATA_HPP_


/**
 * These interfaces have to be supported by the data
 * which is send/recv'd by the simulations
 */
class PararealData
{
public:
	/**
	 * send data to given MPI rank
	 */
	virtual void send(
			int i_mpi_rank,
			int i_mpi_comm
	) = 0;

	/**
	 * receive data from given MPI rank
	 */
	virtual void recv(
			int i_mpi_rank,
			int i_mpi_comm
	) = 0;


	/**
	 * Copy operator which has to be overridden by implementing data container
	 */
	virtual
	const PararealData &operator=(
					const PararealData &i_data
	) = 0;


	/**
	 * Deconstructor
	 */
	~PararealData()
	{
	}
};



template <class T>
class PararealDataInherited	:
		public T
{
};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_DATA_HPP_ */
