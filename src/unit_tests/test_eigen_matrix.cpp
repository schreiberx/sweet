/*
 * test_eigen_matrix.cpp
 *
 *  Created on: 21 Nov 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */

#include <sweet/FatalError.hpp>
#include <iostream>
#include <iomanip>
#include <libmath/EigenMatrix3.hpp>


int main(
		int i_argc,
		char *const i_argv[]
)
{

	std::complex<double> eval[3];

	EigenMatrix3 A(1, 0, 0,
				   0, 0, 1,
			       0, 1, 0);

	A.eigen3real(eval);

	cout<<"Eigenvalues of matrix:"<<endl;
	cout<<eval[0]<<endl;
	cout<<eval[1]<<endl;
	cout<<eval[2]<<endl;
	cout<<"-----------------------"<<endl;
}
