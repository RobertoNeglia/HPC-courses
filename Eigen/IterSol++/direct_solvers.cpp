#include "cg.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/src/Core/util/Constants.h>
#include <functional>
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

using std::cout;
using std::endl;

using namespace LinearAlgebra;

using SpMat = Eigen::SparseMatrix<double>;
using Vect = Eigen::VectorXd;

using create_matrix_t = std::function<void(SpMat &, const int &)>;

void hilbert_matrix(SpMat &A, const int &size) {
  // Create Hilbert matrix
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      A.coeffRef(i, j) = 1.0 / (i + j + 1);
    }
  }
}

void laplace_hilbert_distubed_matrix(SpMat &A, const int &size) {
  int n = size / 50;
  // Create Hilbert matrix
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A.coeffRef(i, j) = 1.0 / (i + j + 1);
    }
  }

  // Complete the matrix as a Laplacian
  for (int i = 0; i < size; i++) {
    A.coeffRef(i, i) += 2.0;
    if (i > 0)
      A.coeffRef(i, i - 1) += -1.0;
    if (i < size - 1)
      A.coeffRef(i, i + 1) += -1.0;
  }
}

int main(int argc, char **argv) {

  // Instantiate Hilbert matrix
  int size = 1000;
  SpMat M(size, size);

  // hilbert_matrix(M, size);
  // laplace_hilbert_distubed_matrix(M, size);

  if (argc != 2) {
    cout << "Error" << endl << "Usage: ./file.o matrix_filename.mtx " << endl;
    return 0;
  }
  std::string matrix_filename(argv[1]);

  Eigen::loadMarket(M, matrix_filename);

  cout << "M - Mt = " << (M - SpMat(M.transpose())).norm() << endl;

  size = M.rows();
  double tol = 1.e-12;
  int result, maxit = 1000;

  // Create right hand size b
  Vect xe = Vect::Ones(size);
  Vect b = M * xe;
  Vect x(size);

  Eigen::DiagonalPreconditioner<double> D(M);

  // Eigen built-in CG method
  Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(M);
  if (cg.info() != Eigen::Success) {
    cout << "Couldn't compute the matrix" << endl;
    return 0;
  }

  x = cg.solve(b);

  cout << "Solve with Eigen CG method" << endl;
  cout << "Number of iterations: " << cg.iterations() << endl;
  cout << "Estimated error: " << cg.error() << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  // Hand-made CG method
  x = 0 * x;

  result = CG(M, x, b, D, maxit, tol);

  cout << "Solve with hand-made CG method" << endl;
  cout << "Number of iterations: " << maxit << endl;
  cout << "Tolerance achieved: " << tol << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  // Eigen direct methods

  // Start with Cholesky Factorization
  Eigen::SimplicialLDLT<SpMat> Cholesky(M);
  // Compute the factorization
  Cholesky.compute(M);
  if (Cholesky.info() != Eigen::Success) {
    cout << "Couldn't factorize the matrix" << endl;
    return 0;
  }
  // Solve the linear system
  x = Cholesky.solve(b);

  cout << "Solve with Cholesky factorization" << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  // LU Factorization
  Eigen::SparseLU<SpMat> LU(M);
  // Compute the factorization
  LU.compute(M);
  if (LU.info() != Eigen::Success) {
    cout << "Couldn't factorize the matrix" << endl;
    return 0;
  }
  // Solve the linear system
  x = LU.solve(b);

  cout << "Solve with LU factorization" << endl;
  cout << "Effective error: " << (x - xe).norm() << endl;
  cout << endl;

  return 0;
}