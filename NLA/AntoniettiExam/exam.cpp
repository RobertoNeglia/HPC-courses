#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h>
#include <cstdlib> // System includes
#include <iostream>

using std::cout;
using std::endl;

#include "IterSolvers/gmres.hpp"

int main(int argc, char **argv) {
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat = Eigen::SparseMatrix<double>;
  using SpVec = Eigen::VectorXd;

  bool verbose = false;

  int n = 1000;
  SpMat A(n, n);
  A.reserve(n * 4);
  for (int i = 0; i < n; i++) {
    A.coeffRef(i, i) = 2.0 * (i + 1);
    if (i > 0)
      A.coeffRef(i, i - 1) = i;
    if (i < n - 1)
      A.coeffRef(i, i + 1) = (i + 1);

    A.coeffRef(i, n - i - 1) += -(i + 1);
  }

  cout << "Matrix size: " << n << "x" << n << endl;
  cout << "Non-zero entries: " << A.nonZeros() << endl;

  if (verbose)
    cout << A << endl;

  SpVec xstar = SpVec::Ones(n);
  //   SpVec xstar = SpVec::Constant(n, 1);
  if (verbose)
    cout << xstar << endl;

  // Creating right hand size b
  SpVec b = A * xstar;

  // Define static parameters
  int maxiter = n;
  int restart = 1000;
  double tol = 1.e-12;

  // Creating the Jacobi (diagonal) preconditioner
  Eigen::DiagonalPreconditioner<double> D(A);

  SpVec x(n);
  // SpVec x = SpVec::Zero(n);

  // Calling GMRES method from gmres.hpp
  int flag = GMRES(A, x, b, D, restart, maxiter, tol);

  cout << "GMRES with no restart return flag: " << flag << endl;
  cout << "Number of iterations: " << maxiter << endl;
  cout << "Tolerance achieved: " << tol << endl;
  cout << "Error: " << (x - xstar).norm() << endl;

  // Calling GMRES method from gmres.hpp with restart = 50
  x = 0 * x;
  maxiter = n;
  tol = 1.e-12;
  restart = 50;
  flag = GMRES(A, x, b, D, restart, maxiter, tol);

  cout << "GMRES with restart = " << restart << " return flag: " << flag
       << endl;
  cout << "Number of iterations: " << maxiter << endl;
  cout << "Tolerance achieved: " << tol << endl;
  cout << "Error: " << (x - xstar).norm() << endl;

  if (verbose)
    cout << x << endl;
  return 0;
}
