#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;
  std::string matrixname = "nos1.mtx";
  SparseMatrix<double> A;
  loadMarket(A, matrixname);

  // Print matrix infos
  cout << "Matrix size: " << A.rows() << "x" << A.cols() << endl;
  cout << "Number of non-zero entries: " << A.nonZeros() << endl;
  SparseMatrix<double> B = SparseMatrix<double>(A.transpose()) - A;
  cout << "Norm of At - A: " << B.norm() << endl;
  if (B.norm() > 1)
    cout << "Not really that symmetric, isn't it" << endl;
  else
    cout << "Matrix is symmetrix" << endl;

  MatrixXd M(A);
  EigenSolver<MatrixXd> esolver(M);
  if (esolver.info() != Eigen::Success)
    abort();
  std::cout << "Eigenvalues of matrix \n" << esolver.eigenvalues() << std::endl;

  SparseMatrix<double> symA = 0.5 * (SparseMatrix<double>(A.transpose()) + A);

  MatrixXd symM(symA);
  SelfAdjointEigenSolver<MatrixXd> saesolver(symM);
  if (saesolver.info() != Eigen::Success)
    abort();
  std::cout << "Eigenvalues of symmetrix part of matrix:\n"
            << saesolver.eigenvalues() << endl;

  return 0;
}