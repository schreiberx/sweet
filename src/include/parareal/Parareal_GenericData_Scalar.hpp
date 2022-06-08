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

	DataContainer<double>* data = nullptr;

public:
	DataContainer<double>* get_pointer_to_data_Scalar() const override
	{
		return this->data;
	};

	void dataArrays_to_GenericData_Scalar(
						double &u
						) override
	{
		this->get_pointer_to_data_Scalar()->simfields[0] = u;
	}

	void GenericData_Scalar_to_dataArrays(
						double &u
						) override
	{
		u = this->get_pointer_to_data_Scalar()->simfields[0];
	}


public:

	Parareal_GenericData_Scalar():
		Parareal_GenericData()
	{
		this->allocate_data();
	}

	Parareal_GenericData_Scalar(Parareal_GenericData_Scalar &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_Scalar());
		for (int i = 0; i < N; i++)
			this->data->simfields[i] = i_data.get_pointer_to_data_Scalar()->simfields[i];
	};

	Parareal_GenericData_Scalar& operator=(const Parareal_GenericData &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_Scalar());
		for (int i = 0; i < N; i++)
			this->data->simfields[i] = i_data.get_pointer_to_data_Scalar()->simfields[i];
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

/////#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	// size in bytes (for MPI)
	// size of each simfield of data
	std::size_t size()
	{
		return 1;
	}

	void serialize(double *data)
	{
		for (int i = 0; i < N; i++)
		{
			///std::memcpy(data + i, &this->data->simfields[i], 1);
			std::copy(&this->data->simfields[i], &this->data->simfields[i + 1], data);
		}
	};

	void deserialize(double *data)
	{
		for (int i = 0; i < N; i++)
		{
			///std::memcpy(&this->data->simfields[i], data + i, 1);
			std::copy(&data[i], &data[i + 1], &this->data->simfields[i]);
		}
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

	double reduce_norm1()
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += std::abs(this->data->simfields[k]);
		return e;
	}

	double reduce_norm2()
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += this->data->simfields[k] * this->data->simfields[k];
		return std::sqrt(e);
	}


	bool check_for_nan()
	{
		bool found_nan = false;
		for (int i = 0; i < N; i++)
			if (std::isnan(this->data->simfields[i]))
			{
				found_nan = true;
				break;
			}
		return found_nan;
	}

	Parareal_GenericData& operator+=(const Parareal_GenericData &i_data)
	{
#if SWEET_PARAREAL
		assert(this->data->time == i_data.get_pointer_to_data_Scalar()->time);
#endif
		assert(this->data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		for (int i = 0; i < N; i++)
			this->data->simfields[i] += i_data.get_pointer_to_data_Scalar()->simfields[i];

		return *this;
	}


	Parareal_GenericData& operator-=(const Parareal_GenericData &i_data)
	{
#if SWEET_PARAREAL
		assert(this->data->time == i_data.get_pointer_to_data_Scalar()->time);
#endif
		assert(this->data->nb_fields == i_data.get_pointer_to_data_Scalar()->nb_fields);

		for (int i = 0; i < N; i++)
			this->data->simfields[i] -= i_data.get_pointer_to_data_Scalar()->simfields[i];

		return *this;
	}

	Parareal_GenericData& operator*=(const double v)
	{

		for (int i = 0; i < N; i++)
			this->data->simfields[i] *= v;

		return *this;
	}


	void physical_print()
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << ": " << this->data->simfields[i] << std::endl;
		}

	}
};

#endif
