/*
 * Parareal_GenericData_Scalar.hpp
 *
 *  Created on: 25 Feb 2022
 *      Authors: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *               Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_SCALAR_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_GENERICDATA_SCALAR_HPP_

#include <assert.h>
#include <parareal/Parareal_GenericData.hpp>


template <int N>
class Parareal_GenericData_Scalar :
		public Parareal_GenericData
{
	class DataContainer_Scalar :
			public Parareal_GenericData::DataContainer<double>
	{

	public:

		DataContainer_Scalar()
		{
			this->simfields = new double[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = 0.;
		};

		DataContainer_Scalar(
				double i_data
		)
		{
			this->simfields = new double[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = i_data;
		};

		DataContainer_Scalar(
				double* i_data
		)
		{
			this->simfields = new double[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = i_data[i];
		};


		DataContainer_Scalar(DataContainer_Scalar &i_data)
		{
			this->time = i_data.time;
			this->level = i_data.level;
			this->nb_fields = i_data.nb_fields;
			for (int i = 0; i < N; i++)
				this->simfields[i] = i_data.simfields[i];
		};
	};


public:

	DataContainer<double>* data;

public:
	DataContainer<double>* get_pointer_to_data_Scalar() const override
	{
		return this->data;
	};


public:

	Parareal_GenericData_Scalar():
		Parareal_GenericData()
	{
	}



//////////	Parareal_GenericData_Scalar(double i_time, int i_level = 0):
//////////		Parareal_GenericData(N, i_time, i_level)
//////////		//data(0)
//////////	{
//////////		this->allocate_data(N, i_time, i_level);
////template <class t_dataType, int N>
//////////	}
//////////
//////////	Parareal_GenericData_Scalar(
//////////			double i_time,
//////////			double i_level,
//////////			double i_data
//////////	)	:
//////////		Parareal_GenericData(i_time, i_level)
//////////		///data(i_data)
//////////	{
//////////		this->allocate_data(N, i_time, i_level);
//////////	}
//////////
//////////	Parareal_GenericData_Scalar(
//////////			double i_time,
//////////			double* i_data,
//////////			double i_level = 0
//////////	)	:
//////////		Parareal_GenericData(i_time, i_level)
//////////	{
//////////		this->allocate_data(N, i_time, i_level);
//////////		///this->data = i_data;
//////////	}

	Parareal_GenericData_Scalar(Parareal_GenericData_Scalar &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_Scalar());
	};

	Parareal_GenericData_Scalar& operator=(const Parareal_GenericData &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_Scalar());
		return *this;
	};

	~Parareal_GenericData_Scalar()
	{
	};

	void allocate_data()
	{
		this->data = new DataContainer_Scalar();
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
		this->allocate_data();
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
					std::abs(this->data->simfields[k]));
		return e;
	}


	Parareal_GenericData& operator+(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_Scalar()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		Parareal_GenericData_Scalar o_data = *this;
		o_data += i_data;

		return o_data;
	}

	void operator+=(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_Scalar()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		for (int i = 0; i < N; i++)
			this->data->simfields[i] += i_data.get_pointer_to_data_Scalar()->simfields[i];

	}

	Parareal_GenericData& operator-(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_Scalar()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		Parareal_GenericData_Scalar o_data = *this;
		o_data -= i_data;

		return o_data;
	}

	void operator-=(const Parareal_GenericData &i_data)
	{
		assert(this->data->time == i_data.get_pointer_to_data_Scalar()->time);
		assert(this->data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		for (int i = 0; i < N; i++)
			this->data->simfields[i] -= i_data.get_pointer_to_data_Scalar()->simfields[i];
	}

/////	Parareal_GenericData<double>& operator*(const Parareal_GenericData<double> &i_data)
/////	{
/////		assert(this->data->time == i_data.data->time);
/////
/////		Parareal_GenericData_Scalar o_data = *this;
/////		o_data *= i_data;
/////
/////		return o_data;
/////	}
/////
/////	void operator*=(const Parareal_GenericData<double> &i_data)
/////	{
/////		assert(this->data->time == i_data.data->time);
/////		assert(this->data->nb_fields = i_data.data->nb_fields);
/////
/////		for (int i = 0; i < N; i++)
/////			this->data->simfields[i] *= i_data.data->simfields[i];
/////	}
/////
/////	Parareal_GenericData<double>& operator/(const Parareal_GenericData<double> &i_data)
/////	{
/////		assert(this->data->time == i_data.data->time);
/////
/////		Parareal_GenericData_Scalar o_data = *this;
/////		o_data /= i_data;
/////
/////		return o_data;
/////	}
/////
/////	void operator/=(const Parareal_GenericData<double> &i_data)
/////	{
/////		assert(this->data->time == i_data.data->time);
/////		assert(this->data->nb_fields = i_data.data->nb_fields);
/////
/////		for (int i = 0; i < N; i++)
/////			this->data->simfields[i] /= i_data.data->simfields[i];
/////	}


};

#endif
