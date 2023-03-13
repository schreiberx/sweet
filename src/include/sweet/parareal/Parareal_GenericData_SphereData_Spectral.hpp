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
#include <sweet/parareal/Parareal_GenericData.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>


namespace sweet
{

template <int N>
class Parareal_GenericData_SphereData_Spectral :
		public Parareal_GenericData
{
	class DataContainer_SphereData_Spectral :
			public Parareal_GenericData::DataContainer<SphereData_Spectral*>
	{

	public:

		DataContainer_SphereData_Spectral(SphereData_Config* i_sphereDataConfig)
		{
			nb_fields = N;
			simfields = new SphereData_Spectral*[N];
			for (int i = 0; i < N; i++)
				simfields[i] = new SphereData_Spectral(i_sphereDataConfig);
		};

		DataContainer_SphereData_Spectral(
				SphereData_Spectral* i_simfields[N]
		)
		{
			nb_fields = N;
			simfields = new SphereData_Spectral*[N];
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_simfields[i]);
		};

		DataContainer_SphereData_Spectral(DataContainer_SphereData_Spectral &i_data)
		{
			nb_fields = N;
			level = i_data.level;
			simfields = new SphereData_Spectral*[N];
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_data.simfields[i]);
		};

		DataContainer_SphereData_Spectral& operator=(const DataContainer_SphereData_Spectral &i_data)
		{
			nb_fields = N;
			level = i_data.level;
			for (int i = 0; i < N; i++)
				*(simfields[i]) = *(i_data.simfields[i]);
			return *this;
		};


		~DataContainer_SphereData_Spectral()
		{
			for (int i = 0; i < N; ++i)
				if (simfields[i])
					delete simfields[i];
			if (simfields)
				delete [] simfields;
		}

	};

public:

	DataContainer<SphereData_Spectral*>* data = nullptr;

public:
	DataContainer<SphereData_Spectral*>* get_pointer_to_data_SphereData_Spectral() const override
	{
		return data;
	};


public:

	Parareal_GenericData_SphereData_Spectral():
		Parareal_GenericData()
	{
		///allocate_data();
	}


	Parareal_GenericData_SphereData_Spectral(Parareal_GenericData_SphereData_Spectral &i_data)
	{
		*(data) = *(i_data.get_pointer_to_data_SphereData_Spectral());
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) = *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]);
		sphereDataConfig = i_data.sphereDataConfig;
	};

	Parareal_GenericData_SphereData_Spectral& operator=(const Parareal_GenericData &i_data)
	override
	{
		*(data) = *(i_data.get_pointer_to_data_SphereData_Spectral());
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) = *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]);
		sphereDataConfig = i_data.sphereDataConfig;
		return *this;
	};


	~Parareal_GenericData_SphereData_Spectral()
	{
		free_data();
	};

	void allocate_data()
	override
	{
		data = new DataContainer_SphereData_Spectral(sphereDataConfig);
	}

	void free_data()
	override
	{
		if (data)
		{
			delete data;
			data = nullptr;
		}
		
	}

#if 0
	/*
	 * M@J: This assignment from double to double* doesn't make sense :-)
	 */
	/**
	 * Setup data
	 */
	void setup(
			double i_data
	)
	{
		data = i_data;
	}
#endif

////#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	// size in bytes (for MPI)
	// size of each simfield of data
	std::size_t size()
	{
		return N * data->simfields[0]->sphereDataConfig->spectral_array_data_number_of_elements;
	}

	void serialize(std::complex<double> *data)
	{
		int s = data->simfields[0]->sphereDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&data->simfields[i]->spectral_space_data[0], &data->simfields[i]->spectral_space_data[s], &data[i * s]);
	};

	void deserialize(std::complex<double> *data)
	{
		int s = data->simfields[0]->sphereDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&data[i * s], &data[(i + 1) * s], &data->simfields[i]->spectral_space_data[0]);
	};

#endif


	void set_time(double i_time)
	override
	{
		data->set_time(i_time);
	}

	double spectral_reduce_maxAbs()
	override
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					data->simfields[k]->spectral_reduce_max_abs());
		return e;
	}

	double spectral_reduce_maxAbs(std::size_t rnorm)
	override
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					data->simfields[k]->spectral_reduce_max_abs(rnorm));
		return e;
	}

	double physical_reduce_maxAbs()
	override
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
					data->simfields[k]->toPhys().physical_reduce_max_abs());
		return e;
	}


	double physical_reduce_norm1()
	override
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += data->simfields[k]->toPhys().physical_reduce_norm1();
		return e;
	}

	double physical_reduce_norm2()
	override
	{
		double e = 0;
		for (int k = 0; k < N; k++)
		{
			double n = data->simfields[k]->toPhys().physical_reduce_norm2();
			e += n * n;
		}
		return std::sqrt(e);
	}

	bool check_for_nan()
	override
	{
		bool found_nan = false;

		int size_n = data->simfields[0]->sphereDataConfig->spectral_modes_n_max;
		int size_m = data->simfields[0]->sphereDataConfig->spectral_modes_m_max;

		for (int i = 0; i < N; i++)
			if (!found_nan)
			{
				for (int m = 0; m < size_m; ++m)
				{
					if (!found_nan)
					{
						for (int n = m; n < size_n; ++n)
							if ( std::isnan(data->simfields[i]->spectral_get_(n, m).real()) ||
								std::isnan(data->simfields[i]->spectral_get_(n, m).imag()) )
							{
								found_nan = true;
								break;
							}
					}
				}
			}

		return found_nan;
	}


	Parareal_GenericData& operator+=(const Parareal_GenericData &i_data)
	override
	{
#if SWEET_PARAREAL
		assert(data->time == i_data.get_pointer_to_data_SphereData_Spectral()->time);
#endif
		assert(data->nb_fields == i_data.get_pointer_to_data_SphereData_Spectral()->nb_fields);
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) += *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]);

		return *this;
	}


	Parareal_GenericData& operator-=(const Parareal_GenericData &i_data)
	override
	{
#if SWEET_PARAREAL
		assert(data->time == i_data.get_pointer_to_data_SphereData_Spectral()->time);
#endif
		assert(data->nb_fields == i_data.get_pointer_to_data_SphereData_Spectral()->nb_fields);
		for (int i = 0; i < N; i++)
			*(data->simfields[i]) -= *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]);

		return *this;
	}

	Parareal_GenericData& operator*=(const double v)
	override
	{

		for (int i = 0; i < N; i++)
			*(data->simfields[i]) *= v;

		return *this;
	}


	void physical_print()
	override
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			data->simfields[i]->toPhys().print();
		}
	}

	void spectral_print()
	override
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			data->simfields[i]->spectral_print();
		}
	}

	void dataArrays_to_GenericData_SphereData_Spectral(
								SphereData_Spectral &phi,
								SphereData_Spectral &vrt,
								SphereData_Spectral &div
							) override
	{
		*(data->simfields[0]) = phi;
		*(data->simfields[1]) = vrt;
		*(data->simfields[2]) = div;
	}

	void GenericData_SphereData_Spectral_to_dataArrays(
								SphereData_Spectral &phi,
								SphereData_Spectral &vrt,
								SphereData_Spectral &div
							) override
	{
		phi = *(data->simfields[0]);
		vrt = *(data->simfields[1]);
		div = *(data->simfields[2]);
	}

	void restrict(const Parareal_GenericData& i_data)
	override
	{
		for (int i = 0; i < N; i++)
			*data->simfields[i] = data->simfields[i]->restrict( *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]) );
	}

	void pad_zeros(const Parareal_GenericData& i_data)
	override
	{
		for (int i = 0; i < N; i++)
			*data->simfields[i] = data->simfields[i]->pad_zeros( *(i_data.get_pointer_to_data_SphereData_Spectral()->simfields[i]) );
	}


};

}

#endif
