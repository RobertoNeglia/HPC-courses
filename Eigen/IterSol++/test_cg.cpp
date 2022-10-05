#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
// #include <cstdlib> // System includes
#include <iostream>

using std::cout;
using std::endl;

#include "cg.hpp"

int main(int argc, char **argv) {
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat = Eigen::SparseMatrix<double>;
  using SpVec = Eigen::VectorXd;

  bool verbose = 1;
  int n = 5;
  SpMat A(n, n); // define matrix
  A.reserve(2998);
  for (int i = 0; i < n; i++) {
    A.coeffRef(i, i) = 2.0 * (i + 1);
    if (i > 0)
      A.coeffRef(i, i - 1) = -i;
    if (i < n - 1)
      A.coeffRef(i, i + 1) = -(i + 1);
  }

  if (verbose)
    cout << A << endl;

  double tol = 1.e-10;      // Convergence tolerance
  int result, maxit = 1000; // Maximum iterations

  std::cout << "Matrix size:" << A.rows() << "X" << A.cols() << std::endl;
  std::cout << "Non zero entries:" << A.nonZeros() << std::endl;

  SpMat B = SpMat(A.transpose()) - A; // Check symmetry
  std::cout << "Norm of A-A.t: " << B.norm() << std::endl;

  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A * e;
  SpVec x(A.rows());
  Eigen::DiagonalPreconditioner<double> D(A); // Create diagonal preconditioner

  // First with eigen CG
  Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(A);
  x = cg.solve(b);
  std::cout << " Eigen native CG" << std::endl;
  std::cout << "#iterations:     " << cg.iterations() << std::endl;
  std::cout << "estimated error: " << cg.error() << std::endl;
  std::cout << "effective error: " << (x - e).norm() << std::endl;

  // Now with hand-made CG
  x = 0 * x;
  result = CG(A, x, b, D, maxit, tol); // Solve system

  std::cout << " hand-made CG " << std::endl;
  cout << "CG flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  std::cout << "Error norm: " << (x - e).norm() << std::endl;

  return result;
}
