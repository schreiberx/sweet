/*
 * PDESWEPlaneMoriZwanzig_NormalModes.hpp
 *
 *  Created on: 13 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_NORMALMODES_HPP_
#define SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_NORMALMODES_HPP_

#include <sweet/core/ErrorBase.hpp>
#include <sweet/expIntegration/ExpFunctions.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <functional>

#if SWEET_EIGEN
#include <Eigen/Eigenvalues>
#endif


#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include "ShackPDESWEPlaneMoriZwanzig.hpp"

/**
 * SWE Plane normal mode
 */
class PDESWEPlaneMoriZwanzigNormalModes
{
public:

	sweet::ErrorBase error;

	typedef double T;
	typedef std::complex<T> complex;

	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackIOData *shackIOData;

	ShackPDESWEPlaneMoriZwanzig *shackPDESWEPlaneMZ;

	double F;
	complex I(0.0, 1.0);

	sweet::ExpFunctions<T> expFunctions;

	bool shackRegistration(sweet::ShackDictionary &io_dict)
	{
		shackPlaneDataOps = io_dict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackTimestepControl = io_dict.getAutoRegistration<sweet::ShackTimestepControl>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackIOData = io_dict.getAutoRegistration<sweet::ShackIOData>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		shackPDESWEPlaneMZ = io_dict.getAutoRegistration<ShackPDESWEPlaneMoriZwanzig>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(io_dict);

		return true;
	}

	bool setup()
	{
		F = shackPDESWEPlaneMZ->F;
		return True;
	}

	void eigendecomposition(
					int i_k0,
					int i_k1,
					complex o_eigenvalues[3],
					complex o_eigenvectors[3][3]
				)
	{
		if (i_k0 == 0 && i_k1 == 0)
		{
			/*
				* http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,0,0%7D,%7B0,0,f%7D,%7B0,-f,0%7D%7D
				* Order need to be changed
				*/
			//R0
			o_eigenvectors[0][0] = 1.;
			o_eigenvectors[1][0] = 0.;
			o_eigenvectors[2][0] = 0.;
			//R-
			o_eigenvectors[0][1] = 0.;
			o_eigenvectors[1][1] = -I.;
			o_eigenvectors[2][1] = 1.;
			//R+
			o_eigenvectors[0][2] = 0.;
			o_eigenvectors[1][2] = I;
			o_eigenvectors[2][2] = 1.;

			o_eigenvalues[1] = 0.;
			o_eigenvalues[1] = -I;
			o_eigenvalues[2] = I;
		}
		else
		{
				/*
					* Compute EV's of
					* Linear operator
					*
					* ///////[     0        -1       bF^{-1/2} ]
					* ///////[     1         0       cF^{-1/2} ]
					* ///////[ bF^{-1/2} cF^{-1/2}       0     ]
					*
					* [     0      bF^{-1/2}  cF^{-1/2} ]
					* [  bF^{-1/2}    0          -1     ]
					* [  cF^{-1/2}    1           0     ]
					*/

			T F = shackPDESWEPlaneMZ->F;
			double FI = 1. / F;
			double Fs = std::sqrt(F);
			double FsI = std::sqrt(FI);
			double normK2 = std::sqrt(i_k0 * i_k0 + i_k1 * i_k1);
			complex b = -i_k0 * I;
			complex c = -i_k1 * I;
			complex om = expFunctions.l_sqrtcplx(FI * normK2 + 1.);

			// R0
			////o_eigenvectors[0][0] = -c * FsI;
			////o_eigenvectors[1][0] = -b * FsI;
			////o_eigenvectors[2][0] = 1.0;
			o_eigenvectors[0][0] = 1. / (b * FsI);
			o_eigenvectors[1][0] = - c / b;
			o_eigenvectors[2][0] = 1.0;

			// R-
			///o_eigenvectors[0][1] = I * (om * i_k0 + i_k1) * FsI / normK2;
			///o_eigenvectors[1][1] = -I * (om * i_k0 * F + i_k1 * (normK2 + F) ) * FsI / (om * normK2);
			///o_eigenvectors[2][1] = 1.0;
			o_eigenvectors[0][1] = I * (om * i_k1 + i_k2) * Fs / ( i_k1 * i_k1 + F ) ;
			o_eigenvectors[1][1] = ( om * i_k0 * i_k1 + normK2 + F ) / ( om * (i_k1 * i_k1 + F) )
			o_eigenvectors[2][1] = 1.0;

			// R+
			///o_eigenvectors[0][2] = -I * (om * i_k0 - i_k1) * FsI / normK2;
			///o_eigenvectors[1][2] = -I * (om * i_k0 * F - i_k1 * (normK2 + F) ) * FsI / (om * normK2);
			///o_eigenvectors[2][2] = 1.0;
			o_eigenvectors[0][2] = I * (-om * i_k1 + i_k2) * Fs / ( i_k1 * i_k1 + F ) ;
			o_eigenvectors[1][2] = ( om * i_k0 * i_k1 - normK2 - F ) / ( om * (i_k1 * i_k1 + F) )
			o_eigenvectors[2][2] = 1.0;

			o_eigenvalues[0] = 0.0;
			o_eigenvalues[1] = -om;
			o_eigenvalues[2] =  om;
		}

		// normalize eigenvectors
		for (int i = 0; i < 3; i++)
		{
			double norm = 0;
			for (int j = 0; j < 3; j++)
			{
				double a = std::abs(o_eigenvectors[j][i]);
				norm += a * a;
			}
			norm = std::sqrt(norm);
			for (int j = 0; j < 3; j++)
				o_eigenvectors[j][i] /= norm;
		}

	}

};

#endif
