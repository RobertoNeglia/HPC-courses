#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>

#include <Eigen/Eigen>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

const auto
timeit(const std::function<void()> &f) {
  using namespace std::chrono;
  const auto start = high_resolution_clock::now();
  f();
  const auto end = high_resolution_clock::now();
  return duration_cast<milliseconds>(end - start).count();
}

void
power_method(const Matrix &A, Vector &b, size_t max_iter, double toll) {
  size_t iter = 0;
  double res  = 1.0;
  Vector v;

    while (iter < max_iter && res > toll) {
      v   = (A * b).normalized();
      res = (v - b).norm();
      b   = v;
      iter++;
    }

  std::cout << "Number of iterations: " << iter << std::endl;
  std::cout << "Tolerance achieved: " << res << std::endl;

  return;
}

double
rayleigh_quotient(const Matrix &A, const Vector &v) {
  double eigenvalue;
  Vector v_t(v.transpose());

  eigenvalue = v_t.dot(A * v) / v_t.dot(v);

  return eigenvalue;
}

int
main(int argc, char **argv) {
  const int    dim      = 1000;
  const size_t max_iter = 5000;
  const double toll     = 1.e-6;
  Matrix       M(dim, dim);

    for (int i = 0; i < dim; i++) {
      M.coeffRef(i, i) = (i + 1) * 2;
    }

  //   std::cout << M << std::endl;

  Vector v = Vector::Ones(dim);
  double max_eigenvalue;

  const auto dt = timeit([&]() {
    power_method(M, v, max_iter, toll);
    max_eigenvalue = rayleigh_quotient(M, v);
  });

  std::cout << "Time elapsed: " << dt << "[ms]" << std::endl;
  //   std::cout << "Dominant eigenvector of M: \n" << v << std::endl;
  std::cout << "Maximum eigenvalue of M: " << max_eigenvalue << std::endl;
  return 0;
}