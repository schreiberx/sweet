/*
 * Parareal_GenericData_PlaneData_Spectral.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_PLANEDATA_SPECTRAL_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_PLANEDATA_SPECTRAL_HPP_

#include <assert.h>
#include <parareal/Parareal_GenericData.hpp>

#define SPLITTED_PLANE_DATA 0

template <int N>
class Parareal_GenericData_PlaneData_Spectral :
		public Parareal_GenericData
{
	class DataContainer_PlaneData_Spectral:
			public Parareal_GenericData::DataContainer<PlaneData*>
	{


	public:

		DataContainer_PlaneData_Spectral()
		{
			this->simfields = new PlaneData*[N];
		};

		DataContainer_PlaneData_Spectral(
#if SPLITTED_PLANE_DATA
				PlaneData_Spectral* i_simfields[N]
#else
				PlaneData* i_simfields[N]
#endif
		)
		{
			this->simfields = new PlaneData*[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = i_simfields[i];
		};

		DataContainer_PlaneData_Spectral(DataContainer_PlaneData_Spectral &i_data)
		{
			this->simfields = new PlaneData*[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = i_data.simfields[i];
		};

		~DataContainer_PlaneData_Spectral()
		{
			if (this->simfields)
				delete [] this->simfields;
		}
	};

public:

	DataContainer<PlaneData*>* data;

public:
	DataContainer<PlaneData*>* get_pointer_to_data_PlaneData_Spectral() const override
	{
		return this->data;
	};

public:

	Parareal_GenericData_PlaneData_Spectral():
#if SPLITTED_PLANE_DATA
		Parareal_GenericData()
#else
		Parareal_GenericData()
#endif
	{
	}


	Parareal_GenericData_PlaneData_Spectral(Parareal_GenericData_PlaneData_Spectral &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_PlaneData_Spectral());
	};

	Parareal_GenericData_PlaneData_Spectral& operator=(const Parareal_GenericData &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_PlaneData_Spectral());
		return *this;
	};


	~Parareal_GenericData_PlaneData_Spectral()
	{
		free_data();
	};


	void set_time(double i_time)
	{
		this->data->set_time(i_time);
	}

	void allocate_data()
	{
		this->data = new DataContainer_PlaneData_Spectral();
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
#if SPLITTED_PLANE_DATA
				PlaneData_Spectral* i_simfields[N]
#else
				PlaneData* i_simfields[N]
#endif
	)
	{
		this->get_pointer_to_data_PlaneData_Spectral()->simfields = i_simfields;
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

	double reduce_maxAbs()
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					this->data->simfields[k]->reduce_maxAbs());
		return e;
	}


	Parareal_GenericData& operator+(const Parareal_GenericData &i_data)
	{

		assert(this->data->time == i_data.get_pointer_to_data_PlaneData_Spectral()->time);
		assert(this->data->nb_fields = i_data.get_pointer_to_data_PlaneData_Spectral()->nb_fields);

		Parareal_GenericData_PlaneData_Spectral<N> o_data = *this;
		o_data += i_data;

		return o_data;
	}

	void operator+=(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_PlaneData_Spectral()->time);
		assert(this->data->nb_fields = i_data.get_pointer_to_data_PlaneData_Spectral()->nb_fields);

		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) += *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]);
	}

	Parareal_GenericData& operator-(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_PlaneData_Spectral()->time);
		assert(this->data->nb_fields = i_data.get_pointer_to_data_PlaneData_Spectral()->nb_fields);

		Parareal_GenericData_PlaneData_Spectral<N> o_data = *this;
		o_data -= i_data;

		return o_data;
	}

	void operator-=(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_PlaneData_Spectral()->time);
		assert(this->data->nb_fields = i_data.get_pointer_to_data_PlaneData_Spectral()->nb_fields);

		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) -= *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]);
	}

//////	Parareal_GenericData<PlaneData*>& operator*(const Parareal_GenericData<PlaneData*> &i_data)
//////	{
//////		assert(this->data->time == i_data.data->time);
//////		assert(this->data->nb_fields = i_data.data->nb_fields);
//////
//////		Parareal_GenericData_PlaneData_Spectral o_data = *this;
//////		o_data *= i_data;
//////
//////		return o_data;
//////	}
//////
//////	void operator*=(const Parareal_GenericData<PlaneData*> &i_data)
//////	{
//////		assert(this->data->time == i_data.data->time);
//////		assert(this->data->nb_fields = i_data.data->nb_fields);
//////
//////		for (int i = 0; i < N; i++)
//////			*(this->data->simfields[i]) *= *(i_data.data->simfields[i]);
//////	}
//////
//////	Parareal_GenericData<PlaneData*>& operator/(const Parareal_GenericData<PlaneData*> &i_data)
//////	{
//////		assert(this->data->time == i_data.data->time);
//////		assert(this->data->nb_fields = i_data.data->nb_fields);
//////
//////		Parareal_GenericData_PlaneData_Spectral o_data = *this;
//////		o_data /= i_data;
//////
//////		return o_data;
//////	}
//////
//////	void operator/=(const Parareal_GenericData<PlaneData*> &i_data)
//////	{
//////		assert(this->data->time == i_data.data->time);
//////		assert(this->data->nb_fields = i_data.data->nb_fields);
//////
//////		for (int i = 0; i < N; i++)
//////			*(this->data->simfields[i]) /= *(i_data.data->simfields[i]);
//////	}


};


#endif
