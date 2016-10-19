/*
 * SPHSetup.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SPHSETUP_HPP_
#define SPHSETUP_HPP_


#include <libmath/shtns_inc.hpp>
#include <fftw3.h>
#include <iostream>
#include <sweet/sweetmath.hpp>


class SphereDataConfig
{
	friend class SphereOperators;
	friend class SphereOperatorsComplex;
	friend class SphereData;
	friend class SphereDataComplex;

private:
	shtns_cfg shtns;

	/**
	 * Number of longitudes
	 */
public:
	int spat_num_lon;

	/**
	 * Number of latitudes
	 */
public:
	int spat_num_lat;

	/**
	 * Number of total longitudes and latitudes
	 */
public:
	int spat_num_elems;


	/**
	 * Number of modes
	 */
public:
	int spec_n_max;
	int spec_m_max;

	/**
	 * Number of total mode variations
	 */
	int spec_num_elems;

	/**
	 * Number of elements for SPH which is based on
	 * complex-valued physical data
	 */
	int cplx_spec_num_elems;


	/**
	 * Array with latitude phi angles
	 *
	 * WARNING: Phi is not the phi from SHTNS
	 */
public:
	double *lat;

	/**
	 * Array with mu = cos(phi) values
	 */
public:
	double *lat_gaussian;

	/**
	 * Array with mu = sin(phi) values
	 */
public:
	double *lat_cogaussian;

public:
	SphereDataConfig()	:
		shtns(nullptr),
		spat_num_lon(-1),
		spat_num_lat(-1),
		spat_num_elems(-1),

		spec_n_max(-1),
		spec_m_max(-1),
		spec_num_elems(-1),
		cplx_spec_num_elems(-1),

		lat(nullptr),
		lat_gaussian(nullptr),
		lat_cogaussian(nullptr)
	{
		refCounter()++;
	}


	std::size_t getArrayIndexByModes(
			int n,
			int m
	)	const
	{
		assert(n >= 0);
		assert(n >= m);

//		return (spec_n_max-im)*im + ((im+1)*im)/2+l;
		return (m*(2*spec_n_max-m+1)>>1)+n;
	}



	std::size_t getArrayIndexByModes_Complex(
			int n,
			int m
	)	const
	{
		assert(n >= 0);
		assert(n >= std::abs(m));

		int idx = n*n+(m+n);
		return idx;
	}


	/**
	 * Return indices with N-variables compactly stored
	 */
	std::size_t getArrayIndexByModes_Complex_NCompact(
			int n,
			int m
	)	const
	{
		assert(n >= 0);
		assert(n >= std::abs(m));

		if (m < 0)
		{
			int minc = spec_m_max+m;
			int row_idx = (((minc-1)*minc)>>1) + minc;
			int idx =  row_idx + (n+m);
			return idx;
		}
		else
		{
			int rel_idx = (m*(2*spec_n_max-m+1)>>1)+n;
			return ((spec_n_max*(spec_n_max+1))>>1)+rel_idx;
		}
	}


private:
	void setup_data()
	{
		spat_num_lat = shtns->nlat;
		spat_num_lon = shtns->nphi;
		spat_num_elems = shtns->nspat;

		spec_n_max = shtns->lmax;
		spec_m_max = shtns->mmax;
		spec_num_elems = shtns->nlm;
		cplx_spec_num_elems = (spec_n_max+1)*(spec_n_max+1);

		if (spec_n_max != spec_m_max)
		{
			std::cerr << "only spec_n_max == spec_m_max currently supported!" << std::endl;
			assert(false);
			exit(1);
		}

		/**
		 * Some safety checks to make sure that we really get what we've asked for
		 */

		/**
		 * TEST: iteration over the modes n,m for real-valued physical space
		 */
		{
			int idx = 0;
			for (int m = 0; m <= spec_m_max; m++)
			{

				int test_idx = getArrayIndexByModes(m,m);

				if (test_idx != idx)
				{
					std::cerr << "IDX TEST NOT SUCCESSFUL (real-valued physical transformation)" << std::endl;
					std::cout << "n=" << m << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
					exit(1);
				}

				for (int n = m; n <= spec_n_max; n++)
				{

					int test_idx2 = getArrayIndexByModes(n,m);
					if (test_idx2 != idx)
					{
						std::cerr << "IDX TEST2 NOT SUCCESSFUL (real-valued physical transformation)" << std::endl;
						std::cout << "n=" << n << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
						exit(1);
					}
					idx++;
				}
			}

			if (idx != spec_num_elems)
			{
				std::cerr << "INTERNAL SPH ERROR (real-valued physical transformation)" << std::endl;
				assert(false);
				exit(1);
			}
		}



		/**
		 * TEST: iteration over the modes n,m for complex-valued physical space
		 *
		 * Note: In SHTNS, the m-coefficients are compactly stored for individual n'l
		 */
		{
			int idx = 0;
			for (int n = 0; n <= spec_n_max; n++)
			{
				int test_idx = getArrayIndexByModes_Complex(n,-n);

				if (test_idx != idx)
				{
					std::cerr << "IDX TEST NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
					std::cout << "n=" << n << ", m=" << -n << "     " << idx << ", " << test_idx << std::endl;
					exit(1);
				}


				for (int m = -n; m <= n; m++)
				{
//					std::cout << std::endl;
//					std::cout << "TESTING inner loop n= " << n << ", m = " << m << std::endl;

					int test_idx2 = getArrayIndexByModes_Complex(n,m);

					if (test_idx2 != idx)
					{
						std::cerr << "IDX TEST2 NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
						std::cout << "n=" << n << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
						exit(1);
					}

					idx++;
				}
			}


			if (idx != cplx_spec_num_elems)
			{
				std::cerr << "INTERNAL SPH ERROR" << std::endl;
				assert(false);
				exit(1);
			}
		}


		/**
		 * TEST: iteration over the modes n,m for complex-valued physical space
		 *
		 * This version tests for the n-coefficients compactly stored for individual m's
		 */
		{
			int idx = 0;
			for (int m = -spec_m_max; m <= spec_m_max; m++)
			{
				int test_idx = getArrayIndexByModes_Complex_NCompact(std::abs(m),m);

				if (test_idx != idx)
				{
					std::cerr << "IDX TEST NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
					std::cout << "n=" << m << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
					exit(1);
				}


				for (int n = std::abs(m); n <= spec_n_max; n++)
				{
					int test_idx2 = getArrayIndexByModes_Complex_NCompact(n,m);

					if (test_idx2 != idx)
					{
						std::cerr << "IDX TEST2 NOT SUCCESSFUL (complex-valued physical transformation)" << std::endl;
						std::cout << "n=" << n << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
						exit(1);
					}

					idx++;
				}
			}
		}

		lat = (double*)fftw_malloc(sizeof(double)*shtns->nlat);

		/*
		 * Colatitude is 0 at the north pole and 180 at the south pole
		 *
		 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
		 */
		for (int i = 0; i < shtns->nlat; i++)
			lat[i] = M_PI*0.5 - ::acos(shtns->ct[i]);

		lat_gaussian = (double*)fftw_malloc(sizeof(double)*shtns->nlat);
		for (int i = 0; i < shtns->nlat; i++)
			lat_gaussian[i] = shtns->ct[i];		/// sin(phi) (SHTNS stores cos(phi))

		lat_cogaussian = (double*)fftw_malloc(sizeof(double)*shtns->nlat);
		for (int i = 0; i < shtns->nlat; i++)
			lat_cogaussian[i] = shtns->st[i];	/// cos(phi) (SHTNS stores sin(phi))
	}


public:
	void setup(
			int mmax,
			int nmax,
			int nphi,
			int nlat
	)
	{
		shtns_verbose(1);			// displays informations during initialization.
		shtns_use_threads(0);		// enable multi-threaded transforms (if supported).

		shtns = shtns_create(
				nmax,
				mmax,
				1,
				(shtns_norm)((int)sht_orthonormal /*| SHT_NO_CS_PHASE*/)
			);

		shtns_set_grid(
				shtns,
				// TODO: replace this with sht_gauss
				(shtns_type)(sht_quick_init | SHT_THETA_CONTIGUOUS),
				//sht_gauss | SHT_THETA_CONTIGUOUS,	// use gaussian grid
				0,
				nlat,		// number of latitude grid points
				nphi		// number of longitude grid points
			);

		setup_data();
	}


	int& refCounter()
	{
		static int ref_counter = 0;
		return ref_counter;
	}


	/**
	 * Setup with given modes.
	 * Spatial resolution will be determined automatically
	 */
	void setupAutoPhysicalSpace(
			int i_mmax,		/// longitude modes
			int i_nmax,		/// latitude modes
			int *o_nphi,	/// physical resolution along longitude
			int *o_nlat		/// physical resolution along latitude
	)
	{
		shtns_verbose(1);			// displays informations during initialization.
		shtns_use_threads(0);		// enable multi-threaded transforms (if supported).


		shtns = shtns_create(
				i_nmax,
				i_mmax,
				1,
				sht_orthonormal
			);

		*o_nphi = 0;
		*o_nlat = 0;

		shtns_set_grid_auto(
				shtns,
				// TODO: replace this with sht_gauss
				(shtns_type)(sht_quick_init | SHT_THETA_CONTIGUOUS),
				//sht_gauss | SHT_THETA_CONTIGUOUS,	// use gaussian grid
				0,
				2,		// use order 2
				o_nlat,
				o_nphi
			);

		setup_data();
	}



	void setupAdditionalModes(
			SphereDataConfig *i_sphConfig,
			int i_additional_modes_longitude,
			int i_additional_modes_latitude
	)
	{
		assert(shtns == nullptr);

		setupAutoPhysicalSpace(
				i_sphConfig->spec_m_max + i_additional_modes_longitude,
				i_sphConfig->spec_n_max + i_additional_modes_latitude,
				&spat_num_lon,
				&spat_num_lat
		);
	}



	void shutdown()
	{
		// check if SPHConfig was initialized
		if (shtns == nullptr)
			return;

		shtns_destroy(shtns);
		shtns = nullptr;

		fftw_free(lat);
		lat = nullptr;

		fftw_free(lat_gaussian);
		lat_gaussian = nullptr;

		fftw_free(lat_cogaussian);
		lat_cogaussian = nullptr;

		assert(refCounter() >= 0);

		if (refCounter() == 0)
		{
			fftw_cleanup_threads();
			fftw_cleanup();
		}
	}



	~SphereDataConfig()
	{
		refCounter()--;
		shutdown();
	}
};



#endif /* SPHSETUP_HPP_ */
