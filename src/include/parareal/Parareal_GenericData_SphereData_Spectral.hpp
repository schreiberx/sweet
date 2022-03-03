/*
 * Parareal_GenericData_SphereData_Spectral.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_SPHEREDATA_SPECTRAL_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_SPHEREDATA_SPECTRAL_HPP_

#include <assert.h>
#include <parareal/Parareal_GenericData.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>

template <int N>
class Parareal_GenericData_SphereData_Spectral :
		public Parareal_GenericData
{
	class DataContainer_SphereData_Spectral :
			public Parareal_GenericData::DataContainer<SphereData_Spectral*>
			//public Parareal_GenericData<SphereData_Spectral*>::DataContainer
	{

	public:

		DataContainer_SphereData_Spectral()
		{
			this->simfields = new SphereData_Spectral*[N];
		};

		DataContainer_SphereData_Spectral(
				SphereData_Spectral* i_simfields[N]
		)
		{
			this->simfields = new SphereData_Spectral*[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = i_simfields[i];
		};

		DataContainer_SphereData_Spectral(DataContainer_SphereData_Spectral &i_data)
		{
			this->simfields = new SphereData_Spectral*[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = i_data.simfields[i];
		};

	};

public:

	DataContainer<SphereData_Spectral*>* data;

public:
	DataContainer<SphereData_Spectral*>* get_pointer_to_data_SphereData_Spectral() const override
	{
		return this->data;
	};


public:

	Parareal_GenericData_SphereData_Spectral():
		Parareal_GenericData()
	{
	}

//////	Parareal_GenericData_SphereData_Spectral(double i_time, int i_level = 0)	:
//////		Parareal_GenericData(i_time, i_level),
//////		data(0)
//////	{
//////	}
//////
//////	Parareal_GenericData_SphereData_Spectral(
//////			double i_time,
//////			double i_level,
//////			double i_data
//////	)	:
//////		Parareal_GenericData(i_time, i_level),
//////		data(i_data)
//////	{
//////	}

	Parareal_GenericData_SphereData_Spectral(Parareal_GenericData_SphereData_Spectral &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_SphereData_Spectral());
	};

	Parareal_GenericData_SphereData_Spectral& operator=(const Parareal_GenericData &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_SphereData_Spectral());
		return *this;
	};


	~Parareal_GenericData_SphereData_Spectral()
	{
		free_data();
	};

	void allocate_data()
	{
		this->data = new DataContainer_SphereData_Spectral();
	}

	void free_data()
	{
		if (this->data)
		{
			delete this->data;
			this->data = nullptr;
		}
		
	}

	/**
	 * Setup data
	 */
	void setup(
			double i_data
	)
	{
		this->data = i_data;
	}

#if SWEET_MPI
	// size in bytes (for MPI)
	std::size_t size()
	{
		return 120398123;
	}

	void serialize(void *data)
	{
	};

	void deserialize(void *data)
	{
	};
#endif


	void set_time(double i_time)
	{
		this->data->set_time(i_time);
	}


	double reduce_maxAbs()
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					this->data->simfields[k]->spectral_reduce_max_abs());
		return e;
	}



	Parareal_GenericData& operator+(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_SphereData_Spectral()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_SphereData_Spectral()->nb_fields);

		Parareal_GenericData_SphereData_Spectral<N> o_data = *this;
		o_data += i_data;

		return o_data;
	}

	void operator+=(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_SphereData_Spectral()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_SphereData_Spectral()->nb_fields);
		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) += *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]);
	}

	Parareal_GenericData& operator-(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_SphereData_Spectral()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_SphereData_Spectral()->nb_fields);

		Parareal_GenericData_SphereData_Spectral<N> o_data = *this;
		o_data -= i_data;

		return o_data;
	}

	void operator-=(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_SphereData_Spectral()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_SphereData_Spectral()->nb_fields);
		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) -= *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]);
	}

/////////	Parareal_GenericData<SphereData_Spectral*, N>& operator*(const Parareal_GenericData<SphereData_Spectral*, N> &i_data)
/////////	{
/////////		assert(this->data.time == i_data->data.time);
/////////
/////////		Parareal_GenericData_SphereData_Spectral o_data = *this;
/////////		o_data *= i_data;
/////////
/////////		return o_data;
/////////	}
/////////
/////////	void operator*=(const Parareal_GenericData<SphereData_Spectral*, N> &i_data)
/////////	{
/////////		assert(this->data.time == i_data->data.time);
/////////		for (int i = 0; i < N; i++)
/////////			this->data.simfield[i] *= i_data.simfield[i];
/////////	}
/////////
/////////	Parareal_GenericData<SphereData_Spectral*, N>& operator/(const Parareal_GenericData<SphereData_Spectral*, N> &i_data)
/////////	{
/////////		assert(this->data.time == i_data->data.time);
/////////
/////////		Parareal_GenericData_SphereData_Spectral o_data = *this;
/////////		o_data /= i_data;
/////////
/////////		return o_data;
/////////	}
/////////
/////////	void operator/=(const Parareal_GenericData<SphereData_Spectral*, N> &i_data)
/////////	{
/////////		assert(this->data.time == i_data->data.time);
/////////		for (int i = 0; i < N; i++)
/////////			this->data.simfield[i] /= i_data.simfield[i];
/////////	}


};


#endif
