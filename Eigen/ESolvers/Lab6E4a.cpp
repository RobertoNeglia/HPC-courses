#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;

int main(int argc, char **argv) {
  using std::cout;
  using std::endl;
  std::string matrixname = "stokes_sym.mtx";
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
    cout << "Matrix is symmetric" << endl;

  MatrixXd M(A);
  SelfAdjointEigenSolver<MatrixXd> saesolver(M);
  if (saesolver.info() != Eigen::Success)
    abort();
  std::cout << "Eigenvalues of matrix:\n" << saesolver.eigenvalues() << endl;

  cout << "Changing the matrix a bit..." << endl;

  for (int i = 0; i < M.rows(); i++) {
    M.coeffRef(i, i) += 2.0;
    if (i > 0)
      M.coeffRef(i, i - 1) += -1.0;
    if (i < M.rows() - 1)
      M.coeffRef(i, i + 1) += -1.0;
  }

  MatrixXd C = MatrixXd(M.transpose()) - M;
  cout << "Norm of Mt - M: " << C.norm() << endl;
  if (C.norm() > 1)
    cout << "Not really that symmetric, isn't it" << endl;
  else
    cout << "Matrix is symmetric" << endl;

  SelfAdjointEigenSolver<MatrixXd> solver(M);
  if (solver.info() != Success)
    abort();
  std::cout << "Eigenvalues of modified matrix:\n"
            << solver.eigenvalues() << endl;

  std::string matrixoutname = "matout.mtx";
  saveMarket(M, matrixoutname);

  return 0;
}