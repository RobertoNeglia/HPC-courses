#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <cstdlib> // System includes
#include <iostream>

using std::cout;
using std::endl;

#include "IterSolvers/bcgstab.hpp"
#include "IterSolvers/cgs.hpp"
#include "IterSolvers/gmres.hpp"

int main(int argc, char **argv) {
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat = Eigen::SparseMatrix<double>;
  using SpVec = Eigen::VectorXd;

  int n = 10;
  SpMat A(n, n);
  SpMat B(n, n);
  SpMat C(n, n); // define matrix
  A.reserve(n * 4);
  B.reserve(n * 3);
  C.reserve(n);
  for (int i = 0; i < n; i++) {
    B.coeffRef(i, i) = 2.0 * (i + 1);
    if (i > 1)
      B.coeffRef(i, i - 1) = i;
    if (i < n - 1)
      B.coeffRef(i, i + 1) = (i + 1);
  }

  cout << B << endl;

  // double tol = 1.e-8;      // Convergence tolerance
  // int result, maxit = 100; // Maximum iterations
  // int restart = 100;       // Restart for gmres

  // std::cout << "Matrix size:" << A.rows() << "X" << A.cols() << std::endl;
  // std::cout << "Non zero entries:" << A.nonZeros() << std::endl;
  // SpMat B = SpMat(A.transpose()) - A;
  // std::cout << "Norm of A-A.t: " << B.norm() << std::endl;
  // SpVec e = SpVec::Ones(A.rows());
  // SpVec b = A * e;
  // SpVec x(A.rows());
  // Eigen::DiagonalPreconditioner<double> D(A); // Create diagonal
  // preconditioner

  // // Solve with CGS method
  // x = 0 * x;
  // result = CGS(A, x, b, D, maxit, tol);
  // cout << "CGS   flag = " << result << endl;
  // cout << "iterations performed: " << maxit << endl;
  // cout << "tolerance achieved  : " << tol << endl;
  // cout << "Error:                " << (x - e).norm() << endl;

  // // Solve with BiCGSTAB method
  // x = 0 * x;
  // maxit = 100;
  // tol = 1.e-8;
  // result = BiCGSTAB(A, x, b, D, maxit, tol);
  // cout << "BiCGSTAB   flag = " << result << endl;
  // cout << "iterations performed: " << maxit << endl;
  // cout << "tolerance achieved  : " << tol << endl;
  // cout << "Error:                " << (x - e).norm() << endl;

  // // Solve with GMRES method
  // x = 0 * x;
  // maxit = 100;
  // tol = 1.e-8;
  // result = GMRES(A, x, b, D, restart, maxit, tol);
  // cout << "GMRES   flag = " << result << endl;
  // cout << "iterations performed: " << maxit << endl;
  // cout << "tolerance achieved  : " << tol << endl;
  // cout << "Error:                " << (x - e).norm() << endl;

  return 0;
}
