#include "cg.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <iostream>

using std::cout;
using std::endl;

using namespace LinearAlgebra;

using SpMat = Eigen::SparseMatrix<double>;
using Vect = Eigen::VectorXd;

int main(int argc, char **argv) {

  // Instantiate Hilbert matrix
  int size = 20, n = 1000;
  SpMat hilbert(n, n);
  hilbert.reserve(3 * n + size * (size - 3) - 2);

  // Create Hilbert matrix
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      hilbert.coeffRef(i, j) = 1.0 / (i + j + 1);
    }
  }

  // Complete the matrix as a Laplacian
  for (int i = 0; i < n; i++) {
    hilbert.coeffRef(i, i) += 2.0;
    if (i > 0)
      hilbert.coeffRef(i, i - 1) += -1.0;
    if (i < n - 1)
      hilbert.coeffRef(i, i + 1) += -1.0;
  }

  double tol = 1.e-12;
  int result, maxit = 2000;

  // Create right hand size b
  Vect xe = Vect::Ones(n);
  Vect b = hilbert * xe;
  Vect x(n);

  Eigen::DiagonalPreconditioner<double> D(hilbert);

  // Eigen built-in CG method
  Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(hilbert);

  x = cg.solve(b);

  cout << "Solve with Eigen CG method" << endl;
  cout << "Number of iterations: " << cg.iterations() << endl;
  cout << "Estimated error: " << cg.error() << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  // Hand-made CG method
  x = 0 * x;

  result = CG(hilbert, x, b, D, maxit, tol);

  cout << "Solve with hand-made CG method" << endl;
  cout << "CG method flag: " << result << endl;
  cout << "Number of iterations: " << maxit << endl;
  cout << "Tolerance achieved: " << tol << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  // Eigen direct methods

  // Start with Cholesky Factorization
  Eigen::SimplicialLDLT<SpMat> Cholesky(hilbert);
  // Compute the factorization
  Cholesky.compute(hilbert);
  // Solve the linear system
  x = Cholesky.solve(b);

  cout << "Solve with Cholesky factorization" << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  // LU Factorization
  Eigen::SparseLU<SpMat> LU(hilbert);
  // Compute the factorization
  LU.compute(hilbert);
  // Solve the linear system
  x = LU.solve(b);

  cout << "Solve with LU factorization" << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  return 0;
}