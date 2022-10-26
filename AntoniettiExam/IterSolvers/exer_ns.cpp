#include <cstdlib>                      // System includes
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using std::endl;
using std::cout;

#include "cgs.hpp"
#include "bcgstab.hpp"
#include "gmres.hpp"

int main(int argc, char** argv)
{
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat=Eigen::SparseMatrix<double>;
  using SpVec=Eigen::VectorXd;
      
  int n = 1000;
  double gam = -0.5;    
  SpMat A(n,n);                      // define matrix
  A.reserve(2997);
  for (int i=0; i<n; i++) {
      A.coeffRef(i, i) = 2.0;
      if(i>1) A.coeffRef(i, i-2) = gam;
      if(i<n-1) A.coeffRef(i, i+1) = 1.0;	
  }

  double tol = 1.e-8;                // Convergence tolerance
  int result, maxit = 100;          // Maximum iterations
  int restart = 100;                 // Restart for gmres

  std::cout<<"Matrix size:"<<A.rows()<<"X"<<A.cols()<<std::endl;
  std::cout<<"Non zero entries:"<<A.nonZeros()<<std::endl;
  SpMat B = SpMat(A.transpose()) - A;
  std::cout<<"Norm of A-A.t: "<<B.norm()<<std::endl;
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());
  Eigen::DiagonalPreconditioner<double> D(A);// Create diagonal preconditioner
  
  // Solve with CGS method
  x=0*x;
  result = CGS(A, x, b, D, maxit, tol);      
  cout << "CGS   flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  cout << "Error:                " << (x-e).norm()<< endl;

  // Solve with BiCGSTAB method
  x=0*x; maxit = 100; tol = 1.e-8;
  result = BiCGSTAB(A, x, b, D, maxit, tol);
  cout << "BiCGSTAB   flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  cout << "Error:                " << (x-e).norm()<< endl;

  // Solve with GMRES method
  x=0*x; maxit = 100; tol = 1.e-8;
  result = GMRES(A, x, b, D, restart, maxit, tol);
  cout << "GMRES   flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  cout << "Error:                " << (x-e).norm()<< endl;

  return 0;
}
