#include "IterSolvers/gmres.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <cstdlib> // System includes
#include <iostream>
int main(int argc, char **argv) {
  using namespace LinearAlgebra;
  using SpMat = Eigen::SparseMatrix<double>;
  using SpVec = Eigen::VectorXd;
  int n = 1000;
  SpMat A(n, n); // define matrix
  for (int i = 0; i < n; i++) {
    A.coeffRef(i, i) = 2.0 * (i + 1);
    A.coeffRef(i, n - i - 1) = -i;
    if (i > 0)
      A.coeffRef(i, i - 1) += i;
    if (i < n - 1)
      A.coeffRef(i, i + 1) += i + 1;
  }
  std::cout << "Matrix size: " << A.rows() << "X" << A.cols() << std::endl;
  std::cout << "Non zero entries: " << A.nonZeros() << std::endl;
  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A * e;
  SpVec x(A.rows());
  //   x = 0 * x;
  // Solve with GMRES method with restart
  double tol = 1.e-12;                        // Convergence tolerance
  int result, maxit = 1000;                   // Maximum iterations
  int restart = 50;                           // Restart gmres
  Eigen::DiagonalPreconditioner<double> D(A); // Create diagonal preconditioner
  result = GMRES(A, x, b, D, restart, maxit, tol);
  std::cout << "GMRES with restart " << std::endl;
  std::cout << "iterations performed: " << maxit << std::endl;
  std::cout << "tolerance achieved : " << tol << std::endl;
  std::cout << "Error: " << (x - e).norm() << std::endl;
  // Solve with GMRES method without restart
  x = 0 * x;
  restart = 1000;
  maxit = 1000;
  tol = 1.e-12;
  result = GMRES(A, x, b, D, restart, maxit, tol);
  std::cout << "GMRES without restart " << std::endl;
  std::cout << "iterations performed: " << maxit << std::endl;
  std::cout << "tolerance achieved : " << tol << std::endl;
  std::cout << "Error norm: " << (x - e).norm() << std::endl;
  return result;
}