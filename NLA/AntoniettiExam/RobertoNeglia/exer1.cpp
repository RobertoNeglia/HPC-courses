#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

#include "cg.hpp"
#include "grad.hpp"

using Matrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;

using std::cout;
using std::endl;

int main() {
  using namespace LinearAlgebra;
  std::string matrix_name("diffreact.mtx");
  Matrix A;
  // Load matrix
  Eigen::loadMarket(A, matrix_name);
  // Check if it's symmetric
  Matrix diff = Matrix(A.transpose()) - A;
  if (diff.norm() < 1.e-6)
    cout << "Matrix is symmetric" << endl;
  else
    cout << "Matrix is not symmetric" << endl;

  cout << "Norm of A: " << A.norm() << endl;
  Matrix Ass = (A - Matrix(A.transpose()));
  cout << "Norm of skew symmetric A: " << Ass.norm() << endl;

  Vector xe = Vector::Ones(A.rows());
  Vector b = A * xe;
  Vector x(A.rows());

  cout << "Norm of b: " << b.norm() << endl;

  Eigen::DiagonalPreconditioner<double> D(A);

  int max_iter = 10000;
  double tol = 1.e-8;

  // Solve Ax = b with grad.hpp
  x = 0 * x;
  int flag = GRAD(A, x, b, D, max_iter, tol);

  cout << "Return flag: " << flag << endl;
  cout << "Number of iterations: " << max_iter << endl;
  cout << "Tolerance achieved: " << tol << endl;
  cout << "Residual error norm: " << (x - xe).norm() << endl;

  // Solve Ax = b with cg.hpp
  max_iter = 10000;
  tol = 1.e-8;
  x = 0 * x;
  flag = CG(A, x, b, D, max_iter, tol);

  cout << "Return flag: " << flag << endl;
  cout << "Number of iterations: " << max_iter << endl;
  cout << "Tolerance achieved: " << tol << endl;
  cout << "Residual error norm: " << (x - xe).norm() << endl;

  return 0;
}