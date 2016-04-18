/*
 * PararealController.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_HPP_


#include <parareal/Parareal_SimulationInstance.hpp>
#include <parareal/Parareal_SimulationVariables.hpp>

#include <iostream>
#include <fstream>
#include <string>


/**
 * This class takes over the control and
 * calls methods offered via PararealSimulation.
 *
 * \param t_SimulationInstance	class which implements the Parareal_SimulationInstance interfaces
 */
template <class t_SimulationInstance>
class Parareal_Controller_Serial
{
public:
	/**
	 * Array with instantiations of PararealSimulations
	 */
	Parareal_SimulationInstance *pararealSimulations = nullptr;

	PararealSimulationVariables *pVars;


	/*
	 * Prefix std::cout with string
	 *
	 * Source: http://stackoverflow.com/questions/27336335/c-cout-with-prefix
	 */
	class prefixbuf
	    : public std::streambuf
	{
	    std::string     prefix;
	    std::streambuf* sbuf;
	    bool            need_prefix;

	    int sync() {
	        return this->sbuf->pubsync();
	    }

	    int overflow(int c) {
	        if (c != std::char_traits<char>::eof()) {
	            if (this->need_prefix
	                && !this->prefix.empty()
	                && this->prefix.size() != this->sbuf->sputn(&this->prefix[0], this->prefix.size())) {
	                return std::char_traits<char>::eof();
	            }
	            this->need_prefix = c == '\n';
	        }
	        return this->sbuf->sputc(c);
	    }

	public:
	    prefixbuf()
	        : prefix(0)
	        , sbuf(0)
	        , need_prefix(true) {
	    }

	public:
	    void setPrefix(std::string const& i_prefix)
	    {
	    	prefix = i_prefix;
	    }
/*
	public:
	    void setPrefix(int const& i_prefix)
	    {
			std::ostringstream ss;
			ss << i_prefix;
	    	prefix = ss.str();
	    }*/

	public:
	    void setStreamBuf(std::streambuf* i_sbuf)
	    {
	    	sbuf = i_sbuf;
	    }
	};

	prefixbuf prefixBuffer;

	std::streambuf *coutbuf;


	void COUT_PrefixStart(int i_number)
	{
		std::ostringstream ss;
		ss << "[" << i_number << "]";
		prefixBuffer.setPrefix(ss.str());

		std::cout.rdbuf(&prefixBuffer);
	}

	void COUT_PrefixEnd()
	{
		std::cout.rdbuf(coutbuf);
	}

	Parareal_Controller_Serial()
	{
		prefixBuffer.setStreamBuf(coutbuf);
		coutbuf = std::cout.rdbuf();
	}

	void setup(
			PararealSimulationVariables *i_pararealSimVars
	)
	{
		pVars = i_pararealSimVars;

		if (!pVars->enabled)
			return;

		if (pVars->coarse_slices <= 0)
		{
			std::cerr << "Invalid number of coarse slices" << std::endl;
			exit(1);
		}

		// allocate simulation instances
		pararealSimulations = new t_SimulationInstance[pVars->coarse_slices];

		/*
		 * SETUP time frame
		 */
		// size of coarse time step
		double coarse_timestep_size = i_pararealSimVars->max_simulation_time / pVars->coarse_slices;

		COUT_PrefixStart(0);
		pararealSimulations[0].sim_set_timeframe(0, i_pararealSimVars->max_simulation_time);

		for (int k = 1; k < pVars->coarse_slices-1; k++)
		{
			COUT_PrefixStart(1);
			pararealSimulations->sim_set_timeframe(coarse_timestep_size*k, coarse_timestep_size*(k+1));
		}

		COUT_PrefixStart(pVars->coarse_slices-1);
		pararealSimulations[pVars->coarse_slices-1].sim_set_timeframe(i_pararealSimVars->max_simulation_time-coarse_timestep_size, i_pararealSimVars->max_simulation_time);


		/*
		 * Setup first simulation instance
		 */
		COUT_PrefixStart(0);
		pararealSimulations[0].sim_setup_initial_data();

		COUT_PrefixEnd();
	}

	void run()
	{
		/**
		 * Initial propagation
		 */
		pararealSimulations[0].run_timestep_coarse();
		for (int i = 1; i < pVars->coarse_slices; i++)
		{
			// use coarse time step output data as initial data of next coarse time step
			pararealSimulations[i].sim_set_data(
					pararealSimulations[i-1].get_data_timestep_coarse()
				);
		}


		/**
		 * We run as much parareal iterations as there are coarse slices
		 */
		for (int k = 0; k < pVars->coarse_slices; k++)
		{
			/**
			 * Fine time stepping
			 */
			for (int i = 0; i < pVars->coarse_slices; i++)
			{
				pararealSimulations[i].run_timestep_fine();
				pararealSimulations[i].compute_difference();
			}
		}


		COUT_PrefixEnd();
	}
};





#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_HPP_ */
