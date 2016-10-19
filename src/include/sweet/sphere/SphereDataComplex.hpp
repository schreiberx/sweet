/*
 * SPHDataComplex.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHDATACOMPLEX_HPP_
#define SPHDATACOMPLEX_HPP_

#include <complex>
#include <functional>
#include <array>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include <sweet/sphere/SphereDataConfig.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/MemBlockAlloc.hpp>



class SphereDataComplex
{
public:
	SphereDataConfig *sphConfig;

public:
	std::complex<double> *data_spat;
	std::complex<double> *data_spec;

	bool data_spec_valid;
	bool data_spat_valid;

public:
	SphereDataComplex(
			SphereDataConfig *i_sphConfig
	)	:
		sphConfig(i_sphConfig),
		data_spat(nullptr),
		data_spec(nullptr)
	{
		assert(i_sphConfig != 0);

		setup(i_sphConfig);
	}


public:
	SphereDataComplex(
			const SphereDataComplex &i_sph_data
	)	:
		sphConfig(i_sph_data.sphConfig),
		data_spat(nullptr),
		data_spec(nullptr)
	{
		setup(i_sph_data.sphConfig);

		operator=(i_sph_data);
	}


public:
	SphereDataComplex(
			const SphereData &i_sph_data
	)	:
		sphConfig(i_sph_data.sphConfig),
		data_spat(nullptr),
		data_spec(nullptr)
	{
		setup(i_sph_data.sphConfig);

		operator=(i_sph_data);
	}


	/**
	 * Run validation checks to make sure that the physical and spectral spaces match in size
	 */
public:
	inline void check(SphereDataConfig *i_sphConfig)	const
	{
		assert(sphConfig->spat_num_lat == i_sphConfig->spat_num_lat);
		assert(sphConfig->spat_num_lon == i_sphConfig->spat_num_lon);

		assert(sphConfig->spec_m_max == i_sphConfig->spec_m_max);
		assert(sphConfig->spec_n_max == i_sphConfig->spec_n_max);
	}


public:
	SphereDataComplex& operator=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		assert(sphConfig == i_sph_data.sphConfig);

		physical_fromSphereData(i_sph_data);
		return *this;
	}

public:
	SphereDataComplex& operator=(
			const SphereDataComplex &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		if (i_sph_data.data_spat_valid)
			memcpy(data_spat, i_sph_data.data_spat, sizeof(cplx)*sphConfig->spat_num_elems);

		if (i_sph_data.data_spec_valid)
			memcpy(data_spec, i_sph_data.data_spec, sizeof(cplx)*sphConfig->cplx_spec_num_elems);

		data_spec_valid = i_sph_data.data_spec_valid;
		data_spat_valid = i_sph_data.data_spat_valid;

		return *this;
	}



public:
	void physical_RealToSphereData(
			SphereData &o_sph_data
	)
	{
		check(o_sph_data.sphConfig);

		request_data_physical();

		for (int i = 0; i < sphConfig->spat_num_elems; i++)
			o_sph_data.data_spat[i] = data_spat[i].real();

		o_sph_data.data_spec_valid = false;
		o_sph_data.data_spat_valid = true;
	}

public:
	void physical_ImagToSphereData(
			SphereData &o_sph_data
	)
	{
		check(o_sph_data.sphConfig);

		request_data_physical();

		for (int i = 0; i < sphConfig->spat_num_elems; i++)
			o_sph_data.data_spat[i] = data_spat[i].imag();

		o_sph_data.data_spec_valid = false;
		o_sph_data.data_spat_valid = true;
	}


public:
	void physical_fromSphereData(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		i_sph_data.request_data_physical();

		for (int i = 0; i < sphConfig->spat_num_elems; i++)
			data_spat[i] = i_sph_data.data_spat[i];

		data_spec_valid = false;
		data_spat_valid = true;
	}


public:
	void request_data_spectral()	const
	{
		if (data_spec_valid)
			return;

		assert(data_spat_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		spat_cplx_to_SH(sphConfig->shtns, data_spat, data_spec);

		SphereDataComplex *this_var = (SphereDataComplex*)this;

		this_var->data_spat_valid = false;
		this_var->data_spec_valid = true;
	}


	void request_data_physical()	const
	{
		if (data_spat_valid)
			return;

		assert(data_spec_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		SH_to_spat_cplx(sphConfig->shtns, data_spec, data_spat);

		SphereDataComplex *this_var = (SphereDataComplex*)this;

		this_var->data_spec_valid = false;
		this_var->data_spat_valid = true;
	}


	SphereDataComplex operator+(
			const SphereDataComplex &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx] + i_sph_data.data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SphereDataComplex& operator+=(
			const SphereDataComplex &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			data_spec[idx] += i_sph_data.data_spec[idx];

		data_spec_valid = true;
		data_spat_valid = false;

		return *this;
	}


	SphereDataComplex& operator-=(
			const SphereDataComplex &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			data_spec[idx] -= i_sph_data.data_spec[idx];

		data_spec_valid = true;
		data_spat_valid = false;

		return *this;
	}



	SphereDataComplex operator-(
			const SphereDataComplex &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx] - i_sph_data.data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SphereDataComplex operator-()
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = -data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SphereDataComplex operator*(
			const SphereDataComplex &i_sph_data
	)	const
	{
		check(i_sph_data.sphConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SphereDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int i = 0; i < sphConfig->cplx_spec_num_elems; i++)
			out_sph_data.data_spat[i] = i_sph_data.data_spat[i]*data_spat[i];

		out_sph_data.data_spec_valid = false;
		out_sph_data.data_spat_valid = true;

		return out_sph_data;
	}


	SphereDataComplex operator*(
			const double i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]*i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	const SphereDataComplex& operator*=(
			const double i_value
	)	const
	{
		request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			data_spec[idx] *= i_value;

		return *this;
	}


	const SphereDataComplex& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			data_spec[idx] *= i_value;

		return *this;
	}


	SphereDataComplex operator/(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]/i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SphereDataComplex operator+(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(*this);

		out_sph_data.data_spec[0] += i_value*std::sqrt(4.0*M_PI);

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}



	SphereDataComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		SphereDataComplex out_sph_data(*this);

		out_sph_data.request_data_spectral();
		out_sph_data.data_spec[0] += i_value*std::sqrt(4.0*M_PI);

		return out_sph_data;
	}



	SphereDataComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

		SphereDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->cplx_spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]*i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}



private:
	void setup(SphereDataConfig *i_sphConfig)
	{
		sphConfig = i_sphConfig;

		data_spec_valid = false;
		data_spat_valid = false;

		data_spat = MemBlockAlloc::alloc<cplx>(sphConfig->spat_num_elems * sizeof(cplx));
		data_spec = MemBlockAlloc::alloc<cplx>(sphConfig->cplx_spec_num_elems * sizeof(cplx));
	}

public:
	~SphereDataComplex()
	{
		MemBlockAlloc::free(data_spat, sphConfig->spat_num_elems * sizeof(cplx));
		MemBlockAlloc::free(data_spec, sphConfig->cplx_spec_num_elems * sizeof(cplx));
	}



	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	inline
	SphereDataComplex spectral_solve_helmholtz(
			const std::complex<double> &i_a,
			const std::complex<double> &i_b,
			double r
	)
	{
		SphereDataComplex out(*this);

		const std::complex<double> a = i_a;
		const std::complex<double> b = i_b/(r*r);

		out.spectral_update_lambda(
			[&](
				int n, int m,
				std::complex<double> &io_data
			)
			{
				io_data /= (a + (-b*(double)n*((double)n+1.0)));
			}
		);

		return out;
	}


	inline
	void spectral_update_lambda(
			std::function<void(int,int,cplx&)> i_lambda
	)
	{
		if (data_spat_valid)
			request_data_spectral();


#pragma omp parallel for
		for (int n = 0; n <= sphConfig->spec_n_max; n++)
		{
			int idx = sphConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				i_lambda(n, m, data_spec[idx]);
				idx++;
			}
		}

		data_spat_valid = false;
		data_spec_valid = true;
	}



	void physical_update_lambda_sinphi_grid(
			std::function< void(double,double,std::complex<double>&) > i_lambda	///< lambda function to return value for lat/mu
	)
	{
		physical_update_lambda_gaussian_grid(i_lambda);
	}

	void physical_update_lambda_cosphi_grid(
			std::function<void(double,double,std::complex<double>&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		physical_update_lambda_cogaussian_grid(i_lambda);
	}



	inline
	const std::complex<double>& spectral_get(
			int in,
			int im
	)	const
	{
		static const std::complex<double> zero = {0,0};

		assert(data_spec_valid);

		if (in < 0)
			return zero;

		if (in > sphConfig->spec_n_max)
			return zero;

		if (std::abs(im) > sphConfig->spec_m_max)
			return zero;

		if (std::abs(im) > in)
			return zero;

		assert (im <= sphConfig->spec_m_max);
		return data_spec[sphConfig->getArrayIndexByModes_Complex(in, im)];
	}



	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */
public:
	void physical_update_lambda(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		if (data_spec_valid)
			request_data_physical();

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphConfig->spat_num_lat; j++)
			{
				//double colatitude = acos(shtns->ct[j]);

				/*
				 * Colatitude is 0 at the north pole and 180 at the south pole
				 *
				 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
				 */
				//double lat_degree = M_PI*0.5 - colatitude;
				double lat_degree = sphConfig->lat[j];

				i_lambda(lon_degree, lat_degree, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}


	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-1;1])
	 */
	void physical_update_lambda_gaussian_grid(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		if (data_spec_valid)
			request_data_physical();

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphConfig->spat_num_lat; j++)
			{
				double mu = sphConfig->lat_gaussian[j];

				i_lambda(lon_degree, mu, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}




	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-1;1])
	 */

	void physical_update_lambda_cogaussian_grid(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		if (data_spec_valid)
			request_data_physical();

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
			double lon_degree = (((double)i)/(double)sphConfig->spat_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphConfig->spat_num_lat; j++)
			{
				double comu = sphConfig->shtns->st[j];
				/*
				 * IDENTITAL FORMULATION
				double mu = shtns->ct[j];
				double comu = sqrt(1.0-mu*mu);
				*/

				i_lambda(lon_degree, comu, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}



	/*
	 * Set all values to zero
	 */
	void physical_set_zero()
	{
#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
			for (int j = 0; j < sphConfig->spat_num_lat; j++)
				data_spat[i*sphConfig->spat_num_lat + j] = 0;

		data_spat_valid = true;
		data_spec_valid = false;
	}


	/**
	 * Return the maximum error norm
	 */
	double physical_reduce_error_max(
			const SphereDataComplex &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		double error = -1;

		for (int j = 0; j < sphConfig->spat_num_elems; j++)
		{
			error = std::max(
						error,
						std::abs(
								data_spat[j].real() - i_sph_data.data_spat[j].real()
							)
						);

			error = std::max(
						error,
						std::abs(
								data_spat[j].imag() - i_sph_data.data_spat[j].imag()
							)
						);
		}
		
		return error;
	}


	/**
	 * Return the maximum error norm
	 */
	double physical_reduce_error_max()
	{
		request_data_physical();

		double error = -1;

		for (int j = 0; j < sphConfig->spat_num_elems; j++)
		{
			error = std::max(
						error,
						std::abs(data_spat[j].real())
						);

			error = std::max(
						error,
						std::abs(data_spat[j].imag())
						);
		}

		return error;
	}

	void spectral_print(
			int i_precision = 8
	)	const
	{
		request_data_spectral();

		std::cout << std::setprecision(i_precision);

		/**
		 * WARNING: This follows a different order contrast to how it is stored
		 */
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				std::size_t idx = sphConfig->getArrayIndexByModes(m, m);
				std::cout << data_spec[idx] << "\t";
			}
			std::cout << std::endl;
		}
	}


	void physical_print(
			int i_precision = 8
	)	const
	{
		request_data_physical();

		std::cout << std::setprecision(i_precision);

#if 0
		for (std::size_t i = 0; i < sphConfig->spat_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			lon_degree = lon_degree/M_PI*180.0;

			std::cout << lon_degree;
			if (i < sphConfig->spat_num_lon-1)
				std::cout << "\t";
		}
		std::cout << std::endl;
#endif

        for (int j = sphConfig->spat_num_lat-1; j >= 0; j--)
        {
#if 0
        		double lat_degree = sphConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		std::cout << lat_degree << "\t";
#endif
        		for (int i = 0; i < sphConfig->spat_num_lon; i++)
        		{
        			std::cout << data_spat[i*sphConfig->spat_num_lat+j];
        			if (i < sphConfig->spat_num_lon-1)
        				std::cout << "\t";
        		}
        		std::cout << std::endl;
        }
	}



public:

	void physical_write_file(
			const char *i_filename,
			const char *i_title = "",
			int i_precision = 8
	)	const
	{
		request_data_physical();

		std::ofstream file(i_filename, std::ios_base::trunc);

		file << std::setprecision(i_precision);
		file << "#TI " << i_title << std::endl;
		file << "#TX Longitude" << std::endl;
		file << "#TY Latitude" << std::endl;

		//file << "lat\\lon\t";
		// Use 0 to make it processable by python
		file << "0\t";

		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
//			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			lon_degree = lon_degree/M_PI*180.0;

			file << lon_degree;
			if (i < sphConfig->spat_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphConfig->spat_num_lat-1; j >= 0; j--)
        {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
        		double lat_degree = sphConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		file << lat_degree << "\t";

        		for (int i = 0; i < sphConfig->spat_num_lon; i++)
        		{
        			file << data_spat[i*sphConfig->spat_num_lat+j];
        			if (i < sphConfig->spat_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}

	void physical_write_file_lon_pi_shifted(
			const char *i_filename,
			std::string i_title = "",
			int i_precision = 8
	)
	{
		request_data_physical();

		std::ofstream file(i_filename, std::ios_base::trunc);

		file << std::setprecision(i_precision);
		file << "#TI " << i_title << std::endl;
		file << "#TX Longitude" << std::endl;
		file << "#TY Latitude" << std::endl;

		//file << "lat\\lon\t";
		// Use 0 to make it processable by python
		file << "0\t";

		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
//			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			lon_degree = (lon_degree-M_PI)/M_PI*180.0;

			file << lon_degree;
			if (i < sphConfig->spat_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphConfig->spat_num_lat-1; j >= 0; j--)
        {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
        		double lat_degree = sphConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		file << lat_degree << "\t";

        		for (int i = 0; i < sphConfig->spat_num_lon; i++)
        		{
        			int ia = i+sphConfig->spat_num_lon/2;
        			if (ia >= sphConfig->spat_num_lon)
        				ia -= sphConfig->spat_num_lon;

        			file << data_spat[ia*sphConfig->spat_num_lat+j];
        			if (i < sphConfig->spat_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}
};




/**
 * operator to support operations such as:
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SphereDataComplex operator*(
		const double i_value,
		const SphereDataComplex &i_array_data
)
{
	return ((SphereDataComplex&)i_array_data)*i_value;
}


/**
 * operator to support operations such as:
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SphereDataComplex operator*(
		const std::complex<double> &i_value,
		const SphereDataComplex &i_array_data
)
{
	return i_array_data*i_value;
}


/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 *
 * Otherwise, we'd have to write it as arrayData-1.5
 *
 */
#if 0
inline
static
SphereDataComplex operator-(
		const double i_value,
		const SphereDataComplex &i_array_data
)
{
	return ((SphereDataComplex&)i_array_data).valueMinusThis(i_value);
//	return -(((SPHDataComplex&)i_array_data).operator-(i_value));
}
#endif
/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */

inline
static
SphereDataComplex operator+(
		const double i_value,
		const SphereDataComplex &i_array_data
)
{
	return ((SphereDataComplex&)i_array_data)+i_value;
}

inline
static
SphereDataComplex operator+(
		const std::complex<double> &i_value,
		const SphereDataComplex &i_array_data
)
{
	return i_array_data+i_value;
}

#endif /* SPHDATA_HPP_ */
