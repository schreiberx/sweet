/*
 * test_eigen_matrix.cpp
 *
 *  Created on: 21 Nov 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */

#include <iostream>
#include <iomanip>
#include <libmath/EigenMatrix3.hpp>
#include <sweet/SWEETError.hpp>

#if SWEET_EIGEN
#include <Eigen/Eigenvalues>
#endif

int main(
		int i_argc,
		char *const i_argv[]
)
{

#if !SWEET_EIGEN
	SWEETError("Cannot test this without Eigen library. Please compile with --eigen=enable");
#endif

	Eigen::MatrixXcf A(3,3) ;
	std::cout<<std::endl;
	std::cout<<"Eigenvalue solver tests"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"----------------------------------------------------------"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Testing real matrix with complex evals 1,-i,i"<<std::endl;

	//Define matrix with ev 1,-i,i
	A(0,0)=1;
	A(0,1)=0;
	A(0,2)=0;
	A(1,0)=0;
	A(1,1)=0;
	A(1,2)=-1;
	A(2,0)=0;
	A(2,1)=1;
	A(2,2)=0;
	std::cout << A << std::endl << std::endl;

	//Eigen solver
	Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
	ces.compute(A);
	//std::cout << "The eigenvalues of A are:" << std::endl << ces.eigenvalues() << std::endl;


	std::complex<double> eval[3];
	std::complex<double> evec[3][3];
	for(int i=0; i<3; i++)
	{
		eval[i]=ces.eigenvalues()[i];
		std::cout << "Eigenvalue "<< i << " : " << eval[i] << std::endl;

	}
	std::cout << "The matrix of eigenvectors, V, is:" << std::endl << ces.eigenvectors() << std::endl << std::endl;

	//std::complex<float> lambda = ces.eigenvalues()[0];
	//std::cout << "Consider the first eigenvalue, lambda = " << lambda << std::endl;
	//Eigen vectors
	//Eigen::VectorXcf v = ces.eigenvectors().col(0);
	//std::cout << "If v is the corresponding eigenvector, then lambda * v = " << std::endl << lambda * v << std::endl;
	//std::cout << "... and A * v = " << std::endl << A * v << std::endl << std::endl;
	//std::cout << "Finally, V * D * V^(-1) = " << std::endl
	//     << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<"----------------------------------------------------------"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Testing complex matrix with complex evals 1,-i,i"<<std::endl;

	std::complex<double> I;;
	I=-1.0;
	I=std::sqrt(I);
	std::cout<<I<<std::endl;
	//Define matrix with ev 1,-i,i
	A(0,0)=1;
	A(0,1)=0;
	A(0,2)=0;
	A(1,0)=0;
	A(1,1)=0;
	A(1,2)=I;
	A(2,0)=0;
	A(2,1)=I;
	A(2,2)=0;
	std::cout << A << std::endl << std::endl;

	//Eigen solver
	ces.compute(A);
	//std::cout << "The eigenvalues of A are:" << std::endl << ces.eigenvalues() << std::endl;

	for(int i=0; i<3; i++)
	{
		eval[i]=ces.eigenvalues()[i];
		std::cout << "Eigenvalue "<< i << " : " << eval[i] << std::endl;

	}
	std::cout << "The matrix of eigenvectors, V, is:" << std::endl << ces.eigenvectors() << std::endl << std::endl;

	//const Eigen::IOFormat fmt(3, Eigen::DontAlignCols, "\t", " ", "", "", "", "");
	//const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", "\n", "", "", "", "");
	//const Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

	//std::complex<float> lambda = ces.eigenvalues()[0];
	//std::cout << "Consider the first eigenvalue, lambda = " << lambda << std::endl;
	//Eigen vectors
	//Eigen::VectorXcf v = ces.eigenvectors().col(0);
	//std::cout << "If v is the corresponding eigenvector, then lambda * v = " << std::endl << lambda * v << std::endl;
	//std::cout << "... and A * v = " << std::endl << A * v << std::endl << std::endl;
	//std::cout << "Finally, V * D * V^(-1) = " << std::endl
	//     << ces.eigenvectors() * ces.eigenvalues().asDiagonal() * ces.eigenvectors().inverse() << std::endl;

	std::cout<< "Done!"<<std::endl;
}
