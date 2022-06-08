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

#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>

class Parareal_GenericData
{

	template <class t_dataType>
	class DataContainer
	{
	public:
		int level;
		double time;
		int nb_fields;

	public:
		t_dataType* simfields;

	public:

		DataContainer()
		{
		};

		DataContainer(int i_nb_fields, double i_time, int i_level = 0) :
			level(i_level),
			time(i_time),
			nb_fields(i_nb_fields)
		{
		};

		DataContainer(DataContainer &i_data) :
			level(i_data.level),
			time(i_data.time),
			nb_fields(i_data.nb_fields)
		{
		};

		DataContainer operator=(const DataContainer &i_data)
		{
			this->level = i_data.level;
			this->time = i_data.time;
			this->nb_fields = i_data.nb_fields;
			return *this;
		};


		~DataContainer()
		{
		}

		void set_time(double i_time)
		{
			this->time = i_time;
		}

	};

public:

	PlaneDataConfig* planeDataConfig = nullptr;
	SphereData_Config* sphereDataConfig = nullptr;

public:
	// different interface functions to avoid template in Parareal_GenericData
	// these interfaces are overridden in the respective child classes
#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
	virtual DataContainer<double>* get_pointer_to_data_Scalar() const
	{
		SWEETError("This interface function should not be called");
		DataContainer<double>* dummy = nullptr;
		return dummy;
	};

	virtual void dataArrays_to_GenericData_Scalar(
							double &u
							)
	{
	};

	virtual void GenericData_Scalar_to_dataArrays(
							double &u
							)
	{
	};


#elif SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
	virtual DataContainer<PlaneData_Spectral*>* get_pointer_to_data_PlaneData_Spectral() const
	{
		SWEETError("This interface function should not be called");
		DataContainer<PlaneData_Spectral*>* dummy = nullptr;
		return dummy;
	};

	virtual void dataArrays_to_GenericData_PlaneData_Spectral(
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
								PlaneData_Spectral &h,
	#endif
								PlaneData_Spectral &u,
								PlaneData_Spectral &v
							)
	{
	};

	virtual void GenericData_PlaneData_Spectral_to_dataArrays(
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
								PlaneData_Spectral &h,
	#endif
								PlaneData_Spectral &u,
								PlaneData_Spectral &v
							)
	{
	};


#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
	virtual DataContainer<SphereData_Spectral*>* get_pointer_to_data_SphereData_Spectral() const
	{
		SWEETError("This interface function should not be called");
		DataContainer<SphereData_Spectral*>* dummy = nullptr;
		return dummy;
	};

	virtual void dataArrays_to_GenericData_SphereData_Spectral(
								SphereData_Spectral &phi,
								SphereData_Spectral &vrt,
								SphereData_Spectral &div
							)
	{
	};

	virtual void GenericData_SphereData_Spectral_to_dataArrays(
								SphereData_Spectral &phi,
								SphereData_Spectral &vrt,
								SphereData_Spectral &div
							)
	{
	};


#endif

	void setup_data_config(PlaneDataConfig* i_planeDataConfig)
	{
		this->planeDataConfig = i_planeDataConfig;
	};

	void setup_data_config(SphereData_Config* i_sphereDataConfig)
	{
		this->sphereDataConfig = i_sphereDataConfig;
	};

public:

	Parareal_GenericData()
	{
	};

	Parareal_GenericData(double i_time, int i_level = 0)
	{
	};

	Parareal_GenericData(Parareal_GenericData &i_data)
	{
	};

	virtual Parareal_GenericData& operator=(const Parareal_GenericData &i_data) = 0;


	//Parareal_GenericData(Parareal_GenericData &&i_data){
	//};

	virtual ~Parareal_GenericData()
	{
	}


	virtual void set_time(double i_time)=0;

	virtual void allocate_data()=0;
	
	virtual void free_data() = 0;

///#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	virtual std::size_t size() = 0;
	////virtual void serialize(void *data) = 0;
	////virtual void deserialize(void *data) = 0;
	#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
	virtual void serialize(double *data) = 0;
	virtual void deserialize(double *data) = 0;
	#else
	virtual void serialize(std::complex<double> *data) = 0;
	virtual void deserialize(std::complex<double> *data) = 0;
	#endif
#endif

	virtual double reduce_maxAbs()=0;
	virtual double reduce_norm1()=0;
	virtual double reduce_norm2()=0;

	virtual bool check_for_nan()=0;

////	virtual Parareal_GenericData operator+(const Parareal_GenericData &i_data) = 0;  // --> not possible because it is an abstract class
	virtual Parareal_GenericData& operator+=(const Parareal_GenericData &i_data) = 0;

////	virtual Parareal_GenericData operator-(const Parareal_GenericData &i_data) = 0;
	virtual Parareal_GenericData& operator-=(const Parareal_GenericData &i_data) = 0;

	virtual Parareal_GenericData& operator*=(const double v) = 0;

	virtual void physical_print() = 0;

};


#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_HPP_ */
