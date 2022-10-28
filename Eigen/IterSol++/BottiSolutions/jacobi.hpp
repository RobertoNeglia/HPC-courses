#include <Eigen/Sparse>
#include <iostream>

namespace LinearAlgebra {

template <class M> void print_matrix(const M &m) {
  using std::cout;
  using std::endl;

  cout << "[";
  for (int i = 0; i < m.rows(); i++) {
    for (int j = 0; j < m.cols(); j++) {
      cout << "\t" << m.coeff(i, j);
    }
    if (i != m.rows() - 1)
      cout << endl;
  }
  cout << "\t]" << endl;
}

template <class matrix, class vector, class preconditioner>
// Jacobi iterative solver implementation
// Function returns 1 if the solution is found before the max number of
// iterations, 0 if not
int jacobi(const matrix &A, vector &x, const vector &b, const preconditioner &p,
           int &maxiter, typename vector::Scalar &tol, bool verbose) {
  using real = typename matrix::Scalar;
  real resid;
  real normb = b.norm();
  vector r = b - A * x;

  if (normb == 0)
    normb = 1;

  resid = r.norm() / normb;

  if (resid < tol) {
    // Already found the solution, end iterations
    maxiter = 0;
    tol = resid;
    return 1;
  }

  for (int i = 0; i < maxiter; i++) {
    x += p.solve(r);
    r = b - A * x;
    resid = r.norm() / normb;
    // if (verbose) {
    //   std::cout << "x" << i << ": " << std::endl;
    //   print_matrix(x);
    //   std::cout << "r" << i << ": " << std::endl;
    //   print_matrix(r);
    // }
    if (resid < tol) {
      maxiter = i;
      tol = resid;
      return 1;
    }
  }

  tol = resid;
  return 0;
}
} // namespace LinearAlgebra