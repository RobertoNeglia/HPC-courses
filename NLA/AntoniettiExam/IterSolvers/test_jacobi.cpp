#include <cstdlib>                      // System includes
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using std::endl;
using std::cout;

#include "cg.hpp"
#include "jacobi.hpp"

int main(int argc, char** argv)
{
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat=Eigen::SparseMatrix<double>;
  using SpVec=Eigen::VectorXd;

  int n = 100;
  SpMat A(n,n);                      // define matrix
  A.reserve(298);
  for (int i=0; i<n; i++) {
      A.coeffRef(i, i) = 2.0;
      if(i>0) A.coeffRef(i, i-1) = -1.0;
      if(i<n-1) A.coeffRef(i, i+1) = -1.0;
  }

  double tol = 1.e-6;                // Convergence tolerance
  int result, maxit = 100;          // Maximum iterations

  std::cout<<"Matrix size:"<<A.rows()<<"X"<<A.cols()<<std::endl;
  std::cout<<"Non zero entries:"<<A.nonZeros()<<std::endl;
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());

  // Solve with CG method
  x=0*x;
  Eigen::DiagonalPreconditioner<double> D(A);// Create diagonal preconditioner
  result = CG(A, x, b, D, maxit, tol);
  cout << "CG   flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  cout << "Error:                " << (x-e).norm()<< endl;

  // Solve with Jacobi method
  x=0*x; maxit = 20000; tol = 1.e-5;
  result = Jacobi(A, x, b, D, maxit, tol);
  cout << "Jacobi flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  cout << "Error:                " << (x-e).norm()<< endl;

  return 0;
}
