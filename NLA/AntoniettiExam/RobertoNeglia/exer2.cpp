#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

using Matrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;

using std::cout;
using std::endl;

int main() {
  int n = 99;
  Matrix A;
  // Load matrix
  for (int i = 0; i < n; i++) {
    A.coeffRef(i, i) = std::abs((n + 1.0) / 2.0 - i + 1) + 1;
    if (i > 0)
      A.coeffRef(i, i - 1) = 0.5;
    if (i < n - 1)
      A.coeffRef(i, i + 1) = 0.5;
  }
  // Check if it's symmetric
  cout << A.coeffRef(0, 0) << endl;
  cout << A.coeffRef(n - 1, n - 1) << endl;
  Eigen::MatrixXd denseA(A);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esolver(A);

  cout << "Eigenvalues of A:" << endl << esolver.eigenvalues() << endl;

  cout << "No time to end the exercise :(" << endl;

  return 0;
}