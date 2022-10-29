#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <Eigen/src/Core/Diagonal.h>
#include <Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h>
#include <Eigen/src/IterativeLinearSolvers/ConjugateGradient.h>
#include <Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h>
#include <iostream>

int main(int argc, char **argv) {
  using matrix = Eigen::SparseMatrix<double>;
  using vector = Eigen::VectorXd;
  using std::cout;
  using std::endl;

  bool verbose = 0;
  int n = 1000;

  // Create TriDiagonal matrix
  matrix A(n, n);
  A.reserve((n * 3) - 2);
  for (int i = 0; i < n; i++) {
    A.coeffRef(i, i) = 2.0 * (i + 1);
    if (i > 0)
      A.coeffRef(i, i - 1) = -i;
    if (i < n - 1)
      A.coeffRef(i, i + 1) = -(i + 1);
  }

  if (verbose)
    cout << A << endl;
  cout << "Matrix size: " << A.rows() << "x" << A.cols() << endl;
  cout << "Number of non-zero entries: " << A.nonZeros() << endl;

  // Symmetry check
  matrix B = Eigen::SparseMatrix<double>(A.transpose()) - A;
  cout << "Norm of A - At = " << B.norm() << endl;

  // CG solver parameters
  double tol = 1.e-10;
  int result_flag, maxit = 1000;

  // Create rhs b
  vector b(n);
  vector t = vector::Ones(n);
  b = A * t;

  vector x(n);

  // Eigen CG calls
  cout << "--- Solving with x = 0 and Diagonal Preconditioner (default)"
       << endl;
  Eigen::ConjugateGradient<matrix, Eigen::Lower | Eigen::Upper> cg;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(A);

  x = cg.solve(b);

  cout << "Number of Iterations: " << cg.iterations() << endl;
  cout << "Estimated error: " << cg.error() << endl;
  cout << "Relative error: " << (x - t).norm() << endl;

  cout << "--- Solving with x = 0 and Incomplete Cholesky Preconditioner"
       << endl;
  Eigen::ConjugateGradient<matrix, Eigen::Lower | Eigen::Upper,
                           Eigen::IncompleteCholesky<double>>
      cgIncChol;
  cgIncChol.setMaxIterations(maxit);
  cgIncChol.setTolerance(tol);
  cgIncChol.compute(A);

  x = cgIncChol.solve(b);

  cout << "Number of Iterations: " << cgIncChol.iterations() << endl;
  cout << "Estimated error: " << cgIncChol.error() << endl;
  cout << "Relative error: " << (x - t).norm() << endl;

  cout << endl;

  cout << "--- Solving with x = 0 and Identity Preconditioner" << endl;
  Eigen::ConjugateGradient<matrix, Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner>
      cgI;
  cgI.setMaxIterations(maxit);
  cgI.setTolerance(tol);
  cgI.compute(A);

  x = cgI.solve(b);

  cout << "Number of Iterations: " << cgI.iterations() << endl;
  cout << "Estimated error: " << cgI.error() << endl;
  cout << "Relative error: " << (x - t).norm() << endl;

  cout << endl;

  cout << "--- Solving with x = 1" << endl;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(A);

  vector x0 = vector::Ones(n);

  x = cg.solveWithGuess(b, x0);

  cout << "Number of Iterations: " << cg.iterations() << endl;
  cout << "Estimated error: " << cg.error() << endl;
  cout << "Relative error: " << (x - t).norm() << endl;

  cout << endl;

  cout << "--- Solving with x = 100 " << endl;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(A);

  vector x1 = vector::Constant(n, 100);

  x = cg.solveWithGuess(b, x1);

  cout << "Number of Iterations: " << cg.iterations() << endl;
  cout << "Estimated error: " << cg.error() << endl;
  cout << "Relative error: " << (x - t).norm() << endl;

  cout << endl;

  cout << "--- Solving with x = 100 and Incomplete Cholesky Preconditioner"
       << endl;
  cgIncChol.setMaxIterations(maxit);
  cgIncChol.setTolerance(tol);
  cgIncChol.compute(A);

  x = cgIncChol.solveWithGuess(b, x1);

  cout << "Number of Iterations: " << cgIncChol.iterations() << endl;
  cout << "Estimated error: " << cgIncChol.error() << endl;
  cout << "Relative error: " << (x - t).norm() << endl;

  cout << endl;

  return 0;
}