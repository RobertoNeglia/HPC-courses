#include "cg.hpp"
#include "cgs.hpp"
#include "mm_readers.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  using namespace Eigen;
  using namespace LinearAlgebra;

  if (argc != 2) {
    std::cerr << " Usage test_load MM_filename" << std::endl;
    return 1;
  }
  std::string matrixFile(argv[1]);
  // Some useful alias
  using SpMat = Eigen::SparseMatrix<double>;
  using SpVec = Eigen::VectorXd;

  // Read matrix
  SpMat A;
  Eigen::loadMarket(A, matrixFile);

  double tol = 1.e-8;                 // Convergence tolerance
  int result, maxit = 1000;           // Maximum iterations for CG
  SpMat B = SpMat(A.transpose()) - A; // Check symmetry
  std::cout << "Norm of A-A.t: " << B.norm() << std::endl;

  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A * e;
  SpVec x(A.rows());
  Eigen::DiagonalPreconditioner<double> D(A);

  // First with Eigen Choleski direct solver
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>
      solver; // LDLT factorization
  solver.compute(A);
  if (solver.info() != Eigen::Success) { // sanity check
    std::cout << "cannot factorize the matrix" << std::endl;
    return 0;
  }

  x = solver.solve(b); // solving
  std::cout << "Solution with Eigen Choleski:" << std::endl;
  std::cout << "effective error: " << (x - e).norm() << std::endl;

  // First with Eigen SparseLU solver
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solvelu; // LU factorization
  solvelu.compute(A);
  if (solvelu.info() != Eigen::Success) { // sanity check
    std::cout << "cannot factorize the matrix"
              << std::endl; // decomposition failed
    return 0;
  }

  x = solvelu.solve(b); // solving
  std::cout << "Solution with Eigen LU:" << std::endl;
  std::cout << "effective error: " << (x - e).norm() << std::endl;

  // Now with hand-made CG
  x = 0 * x;
  result = CG(A, x, b, D, maxit, tol); // Call CG function
  std::cout << "Solution with Conjugate Gradient:" << std::endl;
  std::cout << "iterations performed: " << maxit << std::endl;
  std::cout << "tolerance achieved  : " << tol << std::endl;
  std::cout << "Error norm: " << (x - e).norm() << std::endl;

  return result;
}
