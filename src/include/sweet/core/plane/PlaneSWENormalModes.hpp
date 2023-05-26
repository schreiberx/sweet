/*
 * Staggering.hpp
 *
 *  Created on: 24 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_INCLUDE_SWEET_PLANE_NORMALMODES_HPP_
#define SRC_INCLUDE_SWEET_PLANE_NORMALMODES_HPP_



/**
 * Class which computes eigenvalues and eigenvectors of the linear operator of the SWE on the plane
 * Used by:
 * - Direct exponential integration
 * - Projection operators
 * - etc.
 */
class SWE_Plane_NormalModes
{
public:
	typedef double T;
	typedef std::complex<T> complex;

	sweet::ExpFunctions<T> expFunctions;

	sweet::PlaneData_Config* planeDataConfig = nullptr;

	T s0;
	T s1;

	T f;
	T h;
	T g;

	T sqrt_h;
	T sqrt_g;

public:
	SWE_Plane_NormalModes()
	{
	}

	void setup(
			double i_f,
			double i_h,
			double i_g,
			double i_s0,
			double i_s1,
			sweet::PlaneData_Config* i_planeDataConfig
		)
	{
		this->f = (T)i_f;
		this->h = (T)i_h;
		this->g = (T)i_g;
		this->s0 = (T)i_s0;
		this->s1 = (T)i_s1;

		this->sqrt_h = expFunctions.l_sqrt(this->h);
		this->sqrt_g = expFunctions.l_sqrt(this->g);

		this->planeDataConfig = i_planeDataConfig;
	}

public:
	void eigendecomposition(
					int i_k0,
					int i_k1,
					complex o_eigenvectors[3][3]
				)
	{
		complex dummy1[3];
		complex dummy2[3][3];
		return this->eigendecomposition(
				i_k0, i_k1,
				false,
				false,
				dummy1,
				o_eigenvectors,
				dummy2
			);
	}

public:
	void eigendecomposition(
					int i_k0,
					int i_k1,
					complex o_eigenvectors[3][3],
					complex o_eigenvectors_inv[3][3]
				)
	{
		complex dummy[3];
		return this->eigendecomposition(
				i_k0, i_k1,
				false,
				true,
				dummy,
				o_eigenvectors,
				o_eigenvectors_inv
			);
	}

public:
	void eigendecomposition(
					int i_k0,
					int i_k1,
					complex o_eigenvalues[3],
					complex o_eigenvectors[3][3]
				)
	{
		complex dummy[3][3];
		return this->eigendecomposition(
				i_k0, i_k1,
				true,
				false,
				o_eigenvalues,
				o_eigenvectors,
				dummy
			);
	}

public:
	void eigendecomposition(
					int i_k0,
					int i_k1,
					complex o_eigenvalues[3],
					complex o_eigenvectors[3][3],
					complex o_eigenvectors_inv[3][3]
				)
	{
		return this->eigendecomposition(
				i_k0, i_k1,
				true,
				true,
				o_eigenvalues,
				o_eigenvectors,
				o_eigenvectors_inv
			);
	}

public:
	void eigendecomposition(
					int i_k0,
					int i_k1,
					bool i_compute_eigenvalues,
					bool i_compute_inverse,
					complex o_eigenvalues[3],
					complex o_eigenvectors[3][3],
					complex o_eigenvectors_inv[3][3]
				)
	{

		double k1;
		if (i_k1 < this->planeDataConfig->spectral_data_size[1]/2)
			k1 = (T)i_k1;
		else
			k1 = (T)((int)i_k1-(int)this->planeDataConfig->spectral_data_size[1]);

		double k0 = (T)i_k0;


		complex I(0.0, 1.0);

		complex b = -k0*I;	// d/dx exp(I*k0*x) = I*k0 exp(I*k0*x)
		complex c = -k1*I;

		b = b * this->expFunctions.pi2 / this->s0;
		c = c * this->expFunctions.pi2 / this->s1;

		/*
		 * Matrix with Eigenvectors (column-wise)
		 */
		complex v[3][3];
		complex v_inv[3][3];

		/*
		 * Eigenvalues
		 */
		complex lambda[3];

		if (this->f == 0)
		{
			/*
			 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,0%7D,%7Bg*c,0,0%7D%7D
			 */
			if (i_k0 == 0 && i_k1 == 0)
			{
				v[0][0] = 1;
				v[1][0] = 0;
				v[2][0] = 0;

				v[0][1] = 0;
				v[1][1] = 1;
				v[2][1] = 0;

				v[0][2] = 0;
				v[1][2] = 0;
				v[2][2] = 1;

				if (i_compute_eigenvalues){
					lambda[0] = 0;
					lambda[1] = 0;
					lambda[2] = 0;
				}
			}
			else if (i_k0 == 0)
			{
				v[0][0] = 0;
				v[1][0] = 1;
				v[2][0] = 0;

				v[0][1] = -sqrt_h/sqrt_g;
				v[1][1] = 0;
				v[2][1] = 1;

				v[0][2] = sqrt_h/sqrt_g;
				v[1][2] = 0;
				v[2][2] = 1;

				if (i_compute_eigenvalues){
					lambda[0] = 0;
					lambda[1] = -c * sqrt_g * sqrt_h;
					lambda[2] =  c * sqrt_g * sqrt_h;;
				}
			}
			else if (i_k1 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,0%7D,%7Bg*c*0,0,0%7D%7D
				 */

				v[0][0] = 0;
				v[1][0] = 0;
				v[2][0] = 1;

				v[0][1] = -sqrt_h / sqrt_g;
				v[1][1] = 1;
				v[2][1] = 0;

				v[0][2] = sqrt_h / sqrt_g;
				v[1][2] = 1;
				v[2][2] = 0;

				if (i_compute_eigenvalues){
					lambda[0] = 0;
					lambda[1] = -b * sqrt_g * sqrt_h;
					lambda[2] =  b * sqrt_g * sqrt_h;
				}
			}
			else
			{
				v[0][0] = 0;
				v[1][0] = - c / b;
				v[2][0] = 1.0;

				v[0][1] = -( sqrt_h * expFunctions.l_sqrtcplx( b*b + c*c) ) / (c * sqrt_g);
				v[1][1] = b/c;
				v[2][1] = 1.0;

				v[0][2] =  ( sqrt_h * expFunctions.l_sqrtcplx( b*b + c*c) ) / (c * sqrt_g);
				v[1][2] = b/c;
				v[2][2] = 1.0;

				if (i_compute_eigenvalues){
					lambda[0] = 0.0;
					lambda[1] = -expFunctions.l_sqrtcplx( b*b + c*c) * sqrt_h * sqrt_g;
					lambda[2] =  expFunctions.l_sqrtcplx( b*b + c*c) * sqrt_h * sqrt_g;
				}
			}
		}
		else
		{
			if (i_k0 == 0 && i_k1 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,0,0%7D,%7B0,0,f%7D,%7B0,-f,0%7D%7D
				 */
				v[0][0] = 0;
				v[1][0] = -I;
				v[2][0] = 1;

				v[0][1] = 0;
				v[1][1] = I;
				v[2][1] = 1;

				v[0][2] = 1;
				v[1][2] = 0;
				v[2][2] = 0;

				if (i_compute_eigenvalues){
					lambda[0] =  I * f;
					lambda[1] = -I * f;
					lambda[2] = 0;
				}
			}
			else if (i_k0 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b*0,h*c%7D,%7Bg*b*0,0,f%7D,%7Bg*c,-f,0%7D%7D
				 */
				v[0][0] = f / (c * g);
				v[1][0] = 1;
				v[2][0] = 0;

				v[0][1] = -(c * h) / expFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[1][1] =  -f      / expFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[2][1] = 1;

				v[0][2] = ( c * h) / expFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[1][2] = f        / expFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[2][2] = 1;

				if (i_compute_eigenvalues){
					lambda[0] = 0;
					lambda[1] = -expFunctions.l_sqrtcplx(c*c*g*h - f*f);
					lambda[2] =  expFunctions.l_sqrtcplx(c*c*g*h - f*f);
				}
			}
			else if (i_k1 == 0)
			{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,f%7D,%7Bg*c*0,-f,0%7D%7D
					 */
				v[0][0] = -f / ( b * g);
				v[1][0] = 0;
				v[2][0] = 1;

				v[0][1] = -(b * h) / f;
				v[1][1] =  expFunctions.l_sqrtcplx(-f*f + b*b*g*h) / f;
				v[2][1] = 1;

				v[0][2] = -( b * h) / f;
				v[1][2] = -expFunctions.l_sqrtcplx(-f*f + b*b*g*h) / f;
				v[2][2] = 1;

				if (i_compute_eigenvalues){
					lambda[0] = 0;
					lambda[1] = -expFunctions.l_sqrtcplx(b*b*g*h - f*f);
					lambda[2] =  expFunctions.l_sqrtcplx(b*b*g*h - f*f);
				}
			}
			else
				{
					/*
					 * Compute EV's of
					 * Linear operator
					 *
					 * [ 0  hb  hc ]
					 * [ gb  0   f ]
					 * [ gc -f   0 ]
					 *
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,f%7D,%7Bg*c,-f,0%7D%7D
					 */

					v[0][0] = -f / (b * g);
					v[1][0] = -c / b;
					v[2][0] = 1.0;

					v[0][1] = -( c*f*h + b*h*expFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h) ) /
						( b*c*g*h +    f*expFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h) );
					v[1][1] = -( f*f - b*b*g*h ) / 
						( b*c*g*h +    f*expFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][1] = 1.0;

					v[0][2] = -( -c*f*h + b*h*expFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h) ) / 
						   ( -b*c*g*h + f*expFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h) );
					v[1][2] =  -( -f*f + b*b*g*h ) / 
						   ( -b*c*g*h + f*expFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h) );
					v[2][2] = 1.0;

					/////if (i_k0 == 6 && i_k1 == 11)
					/////{
					/////	std::cout << "AAAAA " << v[0][0] << " " << v[1][0] << " " << v[2][0] << " " << f << " " << b << " " << g << " " << this->expFunctions.pi2 << " " << this->s0 << std::endl;
					/////	std::cout << "AAAAA " << v[0][1] << " " << v[1][1] << " " << v[2][1] << std::endl;
					/////	std::cout << "AAAAA " << v[0][2] << " " << v[1][2] << " " << v[2][2] << std::endl;
					/////}

					if (i_compute_eigenvalues){
						lambda[0] = 0.0;
						lambda[1] = -expFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
						lambda[2] =  expFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
					}
				}
		}


		/*
		 *  Normalize eigenvectors
		 */
		for (int j = 0; j < 3; j++)
		{
			double norm = 0;
			for (int i = 0; i < 3; i++)
				norm += std::abs(v[i][j]) * std::abs(v[i][j]);
			norm = std::sqrt(norm);
			for (int i = 0; i < 3; i++)
				v[i][j] /= norm;
		}

		/*
		 * Invert Eigenvector matrix
		 */

		if (i_compute_inverse){
			v_inv[0][0] =  (v[1][1]*v[2][2] - v[1][2]*v[2][1]);
			v_inv[0][1] = -(v[0][1]*v[2][2] - v[0][2]*v[2][1]);
			v_inv[0][2] =  (v[0][1]*v[1][2] - v[0][2]*v[1][1]);

			v_inv[1][0] = -(v[1][0]*v[2][2] - v[1][2]*v[2][0]);
			v_inv[1][1] =  (v[0][0]*v[2][2] - v[0][2]*v[2][0]);
			v_inv[1][2] = -(v[0][0]*v[1][2] - v[0][2]*v[1][0]);

			v_inv[2][0] =  (v[1][0]*v[2][1] - v[1][1]*v[2][0]);
			v_inv[2][1] = -(v[0][0]*v[2][1] - v[0][1]*v[2][0]);
			v_inv[2][2] =  (v[0][0]*v[1][1] - v[0][1]*v[1][0]);

			complex s = v[0][0]*v_inv[0][0] + v[0][1]*v_inv[1][0] + v[0][2]*v_inv[2][0];

			for (int j = 0; j < 3; j++)
			{
				for (int i = 0; i < 3; i++)
					v_inv[j][i] /= s;
			}

			//Return inverse matrix
			for (int j = 0; j < 3; j++)	{
				for (int i = 0; i < 3; i++)
					o_eigenvectors_inv[j][i] = v_inv[j][i] ;
			}
		}
		if (i_compute_eigenvalues){
			for (int j = 0; j < 3; j++)	{
				o_eigenvalues[j] = lambda[j] ;
			}
		}

		//Return direct matrix
		for (int j = 0; j < 3; j++)	{
			for (int i = 0; i < 3; i++)
				o_eigenvectors[j][i] = v[j][i] ;
		}

		return;

	}


};



#endif
