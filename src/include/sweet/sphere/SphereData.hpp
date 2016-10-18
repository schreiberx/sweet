/*
 * SphereHelper.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHERE_HELPER_HPP_
#define SPHERE_HELPER_HPP_

#include <complex.h>
#include <functional>
#include <array>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <limits>

#include <sweet/sweetmath.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>



class SphereData
{
	friend class SphereDataComplex;

public:
	SphereDataConfig *sphConfig;
//	shtns_cfg shtns;

public:
	std::complex<double> *data_spec;
	double *data_spat;

	bool data_spec_valid;
	bool data_spat_valid;

public:
	SphereData(
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
	SphereData(
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
	SphereData& operator=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		if (i_sph_data.data_spat_valid)
			memcpy(data_spat, i_sph_data.data_spat, sizeof(double)*sphConfig->spat_num_elems);

		if (i_sph_data.data_spec_valid)
			memcpy(data_spec, i_sph_data.data_spec, sizeof(cplx)*sphConfig->spec_num_elems);

		data_spec_valid = i_sph_data.data_spec_valid;
		data_spat_valid = i_sph_data.data_spat_valid;

		return *this;
	}



public:
	const SphereData& spec_copyToDifferentModes(
			SphereData &o_sph_data
	)	const
	{
		request_data_spectral();

		assert (o_sph_data.sphConfig != sphConfig);

		/*
		 *  0 = invalid
		 * -1 = scale down
		 *  1 = scale up
		 */
		int scaling_mode = 0;

		if (sphConfig->spec_m_max < o_sph_data.sphConfig->spec_m_max)
		{
			scaling_mode = 1;
		}
		else if (sphConfig->spec_m_max > o_sph_data.sphConfig->spec_m_max)
		{
			scaling_mode = -1;
		}
//		assert(scaling_mode != 0);


		if (sphConfig->spec_n_max < o_sph_data.sphConfig->spec_n_max)
		{
			assert(scaling_mode != -1);
			scaling_mode = 1;
		}
		else if (sphConfig->spec_n_max > o_sph_data.sphConfig->spec_n_max)
		{
			assert(scaling_mode != 1);
			scaling_mode = -1;
		}
		assert(scaling_mode != 0);


		if (scaling_mode == -1)
		{
			/*
			 * more modes -> less modes
			 */

#pragma omp parallel for
			for (int m = 0; m <= o_sph_data.sphConfig->spec_m_max; m++)
			{
				cplx *dst = &o_sph_data.data_spec[o_sph_data.sphConfig->getArrayIndexByModes(m, m)];
				cplx *src = &data_spec[sphConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(cplx)*(o_sph_data.sphConfig->spec_n_max-m+1);
				memcpy(dst, src, size);
			}
		}
		else
		{
			/*
			 * less modes -> more modes
			 */

			// zero all values
			o_sph_data.spec_set_zero();

#pragma omp parallel for
			for (int m = 0; m <= sphConfig->spec_m_max; m++)
			{
				cplx *dst = &o_sph_data.data_spec[o_sph_data.sphConfig->getArrayIndexByModes(m, m)];
				cplx *src = &data_spec[sphConfig->getArrayIndexByModes(m, m)];

				std::size_t size = sizeof(cplx)*(sphConfig->spec_n_max-m+1);
				memcpy(dst, src, size);
			}
		}

		o_sph_data.data_spat_valid = false;
		o_sph_data.data_spec_valid = true;

		return *this;
	}


	void request_data_spectral()	const
	{
		if (data_spec_valid)
			return;

		assert(data_spat_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		spat_to_SH(sphConfig->shtns, data_spat, data_spec);

		SphereData *this_var = (SphereData*)this;

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
		SH_to_spat(sphConfig->shtns, data_spec, data_spat);

		SphereData *this_var = (SphereData*)this;

		this_var->data_spec_valid = false;
		this_var->data_spat_valid = true;
	}



	SphereData operator+(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx] + i_sph_data.data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}



	SphereData& operator+=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			data_spec[idx] += i_sph_data.data_spec[idx];

		data_spec_valid = true;
		data_spat_valid = false;

		return *this;
	}


	SphereData& operator-=(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			data_spec[idx] -= i_sph_data.data_spec[idx];

		data_spec_valid = true;
		data_spat_valid = false;

		return *this;
	}



	SphereData operator-(
			const SphereData &i_sph_data
	)
	{
		check(i_sph_data.sphConfig);

		request_data_spectral();
		i_sph_data.request_data_spectral();

		SphereData out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx] - i_sph_data.data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SphereData operator-()
	{
		request_data_spectral();

		SphereData out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = -data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SphereData operator*(
			const SphereData &i_sph_data
	)	const
	{
		check(i_sph_data.sphConfig);

		request_data_physical();
		i_sph_data.request_data_physical();

		SphereData out_sph_data(sphConfig);

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_elems; i++)
			out_sph_data.data_spat[i] = i_sph_data.data_spat[i]*data_spat[i];

		out_sph_data.data_spec_valid = false;
		out_sph_data.data_spat_valid = true;

		return out_sph_data;
	}



	SphereData operator*(
			const double i_value
	)	const
	{
		request_data_spectral();

		SphereData out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]*i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}




	const SphereData& operator*=(
			const double i_value
	)	const
	{
		request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			data_spec[idx] *= i_value;

		return *this;
	}





	const SphereData& operator*=(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			data_spec[idx] *= i_value;

		return *this;
	}



	SphereData operator/(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]/i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SphereData operator+(
			double i_value
	)	const
	{
		request_data_spectral();

		SphereData out_sph_data(*this);

		out_sph_data.data_spec[0] += i_value*std::sqrt(4.0*M_PI);

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

		data_spat = MemBlockAlloc::alloc<double>(sphConfig->spat_num_elems * sizeof(double));
		data_spec = MemBlockAlloc::alloc<cplx>(sphConfig->spec_num_elems * sizeof(cplx));
	}


public:
	~SphereData()
	{
		MemBlockAlloc::free(data_spat, sphConfig->spat_num_elems * sizeof(double));
		MemBlockAlloc::free(data_spec, sphConfig->spec_num_elems * sizeof(cplx));
	}


	/**
	 * Solve a Helmholtz problem given by
	 *
	 * (a + b D^2) x = rhs
	 */
	inline
	SphereData spec_solve_helmholtz(
			const double &i_a,
			const double &i_b,
			double r
	)
	{
		SphereData out(*this);

		const double a = i_a;
		const double b = i_b/(r*r);

		out.spec_update_lambda(
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


public:
	/**
	 * Truncate modes which are not representable in spectral space
	 */
	SphereData spat_truncate()
	{
		request_data_physical();

		SphereData out_sph_data(sphConfig);
		spat_to_SH(sphConfig->shtns, data_spat, out_sph_data.data_spec);
		SH_to_spat(sphConfig->shtns, out_sph_data.data_spec, out_sph_data.data_spat);

		out_sph_data.data_spat_valid = true;
		out_sph_data.data_spec_valid = false;

		return out_sph_data;
	}

	/**
	 * Truncate modes which are not representable in spectral space
	 */
	SphereData spec_truncate()
	{
		request_data_spectral();

		SphereData out_sph_data(sphConfig);
		SH_to_spat(sphConfig->shtns, data_spec, out_sph_data.data_spat);
		spat_to_SH(sphConfig->shtns, out_sph_data.data_spat, out_sph_data.data_spec);

		out_sph_data.data_spat_valid = false;
		out_sph_data.data_spec_valid = true;

		return out_sph_data;
	}


	void spec_update_lambda(
			std::function<void(int,int,cplx&)> i_lambda
	)
	{
		if (data_spat_valid)
			request_data_spectral();

#pragma omp parallel for
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			std::size_t idx = sphConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				i_lambda(n, m, data_spec[idx]);
				idx++;
			}
		}

		data_spat_valid = false;
		data_spec_valid = true;
	}



	bool isAnyNaNorInf()
	{
		for (int i = 0; i < sphConfig->spat_num_elems; i++)
		{
			if (std::isnan(data_spat[i]) || std::isinf(data_spat[i]) != 0)
				return true;
		}

		return false;
	}


	const std::complex<double>& spec_get(
			int in,
			int im
	)	const
	{
		static const std::complex<double> zero = {0,0};

		assert(data_spec_valid);

		if (in < 0 ||  im < 0)
			return zero;

		if (in > sphConfig->spec_n_max)
			return zero;

		if (im > sphConfig->spec_m_max)
			return zero;

		if (im > in)
			return zero;

		assert (im <= sphConfig->spec_m_max);
		return data_spec[sphConfig->getArrayIndexByModes(in, im)];
	}


#if 0
	void spec_getElement_enum(
			int n,
			int m,
			std::complex<double> &o_mode_scalar
	)	const
	{
//		if (!data_spec_valid)
//			request_data_spectral();
		assert(data_spec_valid);

		if (	n < 0 ||
				n > sphConfig->spec_n_max	||
				m < 0 ||
				m > sphConfig->spec_m_max	||
				m > n
		)
		{
			o_mode_scalar = {0,0};
			return;
		}

		o_mode_scalar = data_spec[sphConfig->getArrayIndexByModes(n, m)];
	}
#endif


	/*
	 * Set all values to zero
	 */
	void spec_set_zero()
	{
#pragma omp parallel for
		for (int i = 0; i < sphConfig->spec_num_elems; i++)
			data_spec[i] = {0,0};

		data_spat_valid = false;
		data_spec_valid = true;
	}



	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */
	void spat_update_lambda(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
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
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude sin(phi) \in [-1;1])
	 */
	void spat_update_lambda_gaussian_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
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
				double sin_phi = sphConfig->lat_gaussian[j];

				i_lambda(lon_degree, sin_phi, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}




	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters:
	 *   (longitude \in [0;2*pi], Cogaussian latitude cos(phi) \in [0;1])
	 */

	void spat_update_lambda_cogaussian_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
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
				double cos_phi = sphConfig->lat_cogaussian[j];

//				std::cout << cos_phi << std::endl;
				/*
				 * IDENTITAL FORMULATION
				double mu = shtns->ct[j];
				double comu = sqrt(1.0-mu*mu);
				*/

				i_lambda(lon_degree, cos_phi, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}


	void spat_update_lambda_sinphi_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		spat_update_lambda_gaussian_grid(i_lambda);
	}

	void spat_update_lambda_cosphi_grid(
			std::function<void(double,double,double&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		spat_update_lambda_cogaussian_grid(i_lambda);
	}



	/*
	 * Set all values to zero
	 */
	void spat_set_zero()
	{
#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
			for (int j = 0; j < sphConfig->spat_num_lat; j++)
				data_spat[i*sphConfig->spat_num_lat + j] = 0;

		data_spat_valid = true;
		data_spec_valid = false;
	}


	/*
	 * Set all values to a specific value
	 */
	void spat_set_value(
			double i_value
	)
	{
#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
			for (int j = 0; j < sphConfig->spat_num_lat; j++)
				data_spat[i*sphConfig->spat_num_lat + j] = i_value;

		data_spat_valid = true;
		data_spec_valid = false;
	}


	/**
	 * Return the maximum error norm between this and the given data in physical space
	 */
	double spat_reduce_error_max(
			const SphereData &i_sph_data
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
								data_spat[j] - i_sph_data.data_spat[j]
							)
						);
		}
		return error;
	}


	/**
	 * Return the maximum absolute value
	 */
	double spat_reduce_abs_max()
	{
		request_data_physical();

		double error = -1;

		for (int j = 0; j < sphConfig->spat_num_elems; j++)
		{
			error = std::max(
						error,
						std::abs(data_spat[j])
						);
		}
		return error;
	}


	/**
	 * Return the maximum value
	 */
	double spat_reduce_min()
	{
		request_data_physical();

		double error = std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphConfig->spat_num_elems; j++)
		{
			error = std::min(error, data_spat[j]);
		}
		return error;
	}


	/**
	 * Return the minimum value
	 */
	double spat_reduce_max()
	{
		request_data_physical();

		double error = -std::numeric_limits<double>::infinity();

		for (int j = 0; j < sphConfig->spat_num_elems; j++)
		{
			error = std::max(error, data_spat[j]);
		}
		return error;
	}


	void spec_print(int i_precision = 8)	const
	{
		request_data_spectral();

		std::cout << std::setprecision(i_precision);

		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			std::size_t idx = sphConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				std::cout << data_spec[idx] << "\t";
				idx++;
			}
			std::cout << std::endl;
		}
	}


	void spat_print(int i_precision = 8)	const
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


	void spat_write_file(
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

		for (std::size_t i = 0; i < sphConfig->spat_num_lon; i++)
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

	void spat_write_file_lon_pi_shifted(
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

		for (std::size_t i = 0; i < sphConfig->spat_num_lon; i++)
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
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SphereData operator*(
		const double i_value,
		const SphereData &i_array_data
)
{
	return ((SphereData&)i_array_data)*i_value;
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
SphereData operator-(
		const double i_value,
		const SphereData &i_array_data
)
{
	return ((SphereData&)i_array_data).valueMinusThis(i_value);
//	return -(((SPHData&)i_array_data).operator-(i_value));
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

#if 0
inline
static
SphereData operator+(
		const double i_value,
		const SphereData &i_array_data
)
{
	i_array_data.checkConsistency();
	return ((SphereData&)i_array_data)+i_value;
}
#endif

#endif /* SPHDATA_HPP_ */
