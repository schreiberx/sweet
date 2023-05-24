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
#include <sweet/parareal/Parareal_GenericData.hpp>

template <int N>
class Parareal_GenericData_PlaneData_Spectral :
		public sweet::Parareal_GenericData
{
	class DataContainer_PlaneData_Spectral:
			public sweet::Parareal_GenericData::DataContainer<sweet::PlaneData_Spectral*>
	{


	public:

		DataContainer_PlaneData_Spectral(sweet::PlaneData_Config* i_planeDataConfig)
		{
			this->nb_fields = N;
			this->simfields = new sweet::PlaneData_Spectral*[N];
			for (int i = 0; i < N; i++)
				this->simfields[i] = new sweet::PlaneData_Spectral(i_planeDataConfig);
		};

		DataContainer_PlaneData_Spectral(
				sweet::PlaneData_Spectral* i_simfields[N]
		)
		{
			this->nb_fields = N;
			this->simfields = new sweet::PlaneData_Spectral*[N];
			for (int i = 0; i < N; i++)
				*(this->simfields[i]) = *(i_simfields[i]);
		};

		DataContainer_PlaneData_Spectral(DataContainer_PlaneData_Spectral &i_data)
		{
			this->nb_fields = N;
			this->level = i_data.level;
			this->simfields = new sweet::PlaneData_Spectral*[N];
			for (int i = 0; i < N; i++)
				*(this->simfields[i]) = *(i_data.simfields[i]);
		};

		DataContainer_PlaneData_Spectral& operator=(const DataContainer_PlaneData_Spectral &i_data)
		{
			this->nb_fields = N;
			this->level = i_data.level;
			for (int i = 0; i < N; i++)
				*(this->simfields[i]) = *(i_data.simfields[i]);
			return *this;
		};


		~DataContainer_PlaneData_Spectral()
		{
			for (int i = 0; i < N; ++i)
				if (this->simfields[i])
					delete this->simfields[i];
			if (this->simfields)
				delete [] this->simfields;
		}
	};

public:

	DataContainer<sweet::PlaneData_Spectral*>* data = nullptr;

public:
	DataContainer<sweet::PlaneData_Spectral*>* get_pointer_to_data_PlaneData_Spectral() const override
	{
		return this->data;
	};

public:

	Parareal_GenericData_PlaneData_Spectral():
		sweet::Parareal_GenericData()
	{
		////this->allocate_data();
	}

	Parareal_GenericData_PlaneData_Spectral(Parareal_GenericData_PlaneData_Spectral &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_PlaneData_Spectral());
		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]);
		this->planeDataConfig = i_data.planeDataConfig;
	};

	Parareal_GenericData_PlaneData_Spectral& operator=(const Parareal_GenericData &i_data)
	{
		*(this->data) = *(i_data.get_pointer_to_data_PlaneData_Spectral());
		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) = *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]);
		this->planeDataConfig = i_data.planeDataConfig;
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
		this->data = new DataContainer_PlaneData_Spectral(this->planeDataConfig);
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
				sweet::PlaneData_Spectral* i_simfields[N]
	)
	{
		this->get_pointer_to_data_PlaneData_Spectral()->simfields = i_simfields;
	}


////#if SWEET_MPI
#if SWEET_PARAREAL==2 || SWEET_XBRAID
	// size in bytes (for MPI)
	// size of each simfield of data
	std::size_t size()
	{
		return N * this->data->simfields[0]->planeDataConfig->spectral_array_data_number_of_elements;
	}

	void serialize(std::complex<double> *i_data)
	{
		int s = this->data->simfields[0]->planeDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&this->data->simfields[i]->spectral_space_data[0], &this->data->simfields[i]->spectral_space_data[s], &i_data[i * s]);
	};

	void deserialize(std::complex<double> *i_data)
	{
		int s = this->data->simfields[0]->planeDataConfig->spectral_array_data_number_of_elements;
		for (int i = 0; i < N; i++)
			std::copy(&i_data[i * s], &i_data[(i + 1) * s], &this->data->simfields[i]->spectral_space_data[0]);
	};

#endif

	double spectral_reduce_maxAbs()
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
				this->data->simfields[k]->spectral_reduce_max_abs());
		return e;
	}

	double spectral_reduce_maxAbs(std::size_t rnorm)
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
				this->data->simfields[k]->spectral_reduce_max_abs(rnorm));
		return e;
	}

	double physical_reduce_maxAbs()
	{
		double e = -1;
		for (int k = 0; k < N; k++)
			e = std::max( e,
				this->data->simfields[k]->toPhys().physical_reduce_max_abs());
		return e;
	}

	double physical_reduce_norm1()
	{
		double e = 0;
		for (int k = 0; k < N; k++)
			e += this->data->simfields[k]->toPhys().physical_reduce_norm1();
		return e;
	}

	double physical_reduce_norm2()
	{
		double e = 0;
		for (int k = 0; k < N; k++)
		{
			double n = this->data->simfields[k]->toPhys().physical_reduce_norm2();
			e += n * n;
		}
		return std::sqrt(e);
	}


	bool check_for_nan()
	{
		bool found_nan = false;

		int physical_size_x = this->data->simfields[0]->planeDataConfig->physical_data_size[0];
		int physical_size_y = this->data->simfields[0]->planeDataConfig->physical_data_size[1];

		for (int i = 0; i < N; i++)
		{
			sweet::PlaneData_Physical data_phys = this->data->simfields[i]->toPhys();
			if (!found_nan)
			{
				for (int ix = 0; ix < physical_size_x; ++ix)
				{
					if (!found_nan)
					{
						for (int iy = 0; iy < physical_size_y; ++iy)
							if ( std::isnan(data_phys.physical_get(ix, iy)))
							{
								found_nan = true;
								break;
							}
					}
				}
			}
		}

		return found_nan;
	}


	Parareal_GenericData& operator+=(const Parareal_GenericData &i_data)
	{
#if SWEET_PARAREAL
		assert(this->data->time == i_data.get_pointer_to_data_PlaneData_Spectral()->time);
#endif
		assert(this->data->nb_fields == i_data.get_pointer_to_data_PlaneData_Spectral()->nb_fields);

		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) += *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]);

		return *this;
	}

	Parareal_GenericData& operator-=(const Parareal_GenericData &i_data)
	{
#if SWEET_PARAREAL
		assert(this->data->time == i_data.get_pointer_to_data_PlaneData_Spectral()->time);
#endif
		assert(this->data->nb_fields == i_data.get_pointer_to_data_PlaneData_Spectral()->nb_fields);

		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) -= *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]);

		return *this;
	}

	Parareal_GenericData& operator*=(const double v)
	{

		for (int i = 0; i < N; i++)
			*(this->data->simfields[i]) *= v;

		return *this;
	}

	void physical_print()
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			this->data->simfields[i]->toPhys().print();
		}
	}

	void spectral_print()
	{
		for (int i = 0; i < N; i++)
		{
			std::cout << "Field #" << i << std::endl;
			this->data->simfields[i]->spectral_print();
		}
	}


	void dataArrays_to_GenericData_PlaneData_Spectral(
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
								sweet::PlaneData_Spectral &h,
	#endif
								sweet::PlaneData_Spectral &u,
								sweet::PlaneData_Spectral &v
							) override
	{
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
		*(this->data->simfields[0]) = h;
		*(this->data->simfields[1]) = u;
		*(this->data->simfields[2]) = v;
	#elif SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
		*(this->data->simfields[0]) = u;
		*(this->data->simfields[1]) = v;
	#endif
	}

	void GenericData_PlaneData_Spectral_to_dataArrays(
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
								sweet::PlaneData_Spectral &h,
	#endif
								sweet::PlaneData_Spectral &u,
								sweet::PlaneData_Spectral &v
							) override
	{
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
		h = *(this->data->simfields[0]);
		u = *(this->data->simfields[1]);
		v = *(this->data->simfields[2]);
	#elif SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
		u = *(this->data->simfields[0]);
		v = *(this->data->simfields[1]);
	#endif
	}


	void restrict(const Parareal_GenericData& i_data)
	{
		for (int i = 0; i < N; i++)
			*this->data->simfields[i] = this->data->simfields[i]->restrict( *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]) );
			///this->data->simfields[i]->restrict( *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]) );
	}

	void pad_zeros(const Parareal_GenericData& i_data)
	{
		for (int i = 0; i < N; i++)
			*this->data->simfields[i] = this->data->simfields[i]->pad_zeros( *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]) );
			///this->data->simfields[i]->pad_zeros( *(i_data.get_pointer_to_data_PlaneData_Spectral()->simfields[i]) );
	}

};


#endif
