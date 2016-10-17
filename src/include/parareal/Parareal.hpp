/*
 * Parareal.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_HPP_




#if SWEET_PARAREAL==0


/*
 * Create empty Parareal implementations?
 */
class PararealSimulation_Base{};
class PararealData{};
template <class T> class PararealDataInherited	: public T {};

#elif SWEET_PARAREAL==1

#	include <parareal/Parareal_SimulationInstance.hpp>
#	include <parareal/Parareal_Controller_Serial.hpp>
#	include <parareal/Parareal_Data_PlaneData.hpp>


#elif SWEET_PARAREAL==2

#	error "Parareal with MPI not yet supported"

#endif



#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_HPP_ */
