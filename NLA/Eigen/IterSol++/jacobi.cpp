#include <Eigen/Sparse>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <time.h>
#include <unsupported_eigen/Eigen/SparseExtra>

#include "jacobi.hpp"

using std::cout;
using std::endl;

using Matrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;

using create_method_t = std::function<void(Matrix &, const int &)>;

using namespace LinearAlgebra;

void random_matrix(Matrix &m, int size) {
  // Random number seed generator
  srand(time(NULL));

  m.reserve(2 * size);
  int i, j;
  for (int k = 0; k < size; k++) {
    i = rand() % size;
    j = rand() % size;
    m.coeffRef(i, j) = k * 2.32;
    m.coeffRef(k, k) = (k + 4) * 4.67;
  }
}

void tridiagonal_matrix(Matrix &m, int size) {
  srand(time(NULL));

  m.reserve(3 * size);
  int random;

  for (int i = 0; i < size; i++) {
    random = (rand() % 10) - 5;
    if (random != 0) {
      // cout << random << " - ";
      m.coeffRef(i, i) = random * (i + 1) * 3.7;
      if (i > 0)
        m.coeffRef(i, i - 1) = random * 1.2;
      if (i < size - 1)
        m.coeffRef(i, i + 1) = random * 0.6;
    } else
      i--;
  }
}

void load_matrix(Matrix &m, int size) { Eigen::loadMarket(m, "bcsstm12.mtx"); }

void create_matrix(Matrix &m, int size, create_method_t create) {
  create(m, size);
}

int main(int argc, char **argv) {
  bool verbose = 0;

  // Create the Matrix
  int n = 10000; // size of the Matrix
  Matrix A(n, n);

  cout << "Creating matrix..." << endl;
  create_matrix(A, n, tridiagonal_matrix);
  cout << "Matrix created!" << endl;

  n = A.rows();

  if (verbose) {
    print_matrix(A);
  }

  // Right hand side
  Vector b(n);
  // Approximate iteration result Vector
  // Vector xk = Vector::Zero(n);
  Vector xk = Vector::Constant(n, 0);
  // Real solution Vector
  Vector x = Vector::Ones(n);

  // Create right hand side such that x = 1
  b = A * x;

  if (verbose) {
    cout << "xk: " << endl;
    print_matrix(xk);
    cout << "x: " << endl;
    print_matrix(x);
    cout << "b: " << endl;
    print_matrix(b);
  }

  // Create diagonal preconditioner to have D^-1
  Eigen::DiagonalPreconditioner<double> D(A);

  // Jacobi method parameters
  int maxiter = 2000;
  double tol = 1.e-13;

  int result;

  cout << "--------------------------------" << endl;
  {
    using namespace std::chrono;
    auto start = std::chrono::high_resolution_clock::now();
    result = jacobi(A, xk, b, D, maxiter, tol, verbose);
    auto end = std::chrono::high_resolution_clock::now();
    auto dt = duration_cast<milliseconds>(end - start).count();
    cout << "Jacobi method execution time: " << dt << endl;
  }
  cout << "Jacobi result flag: " << result << endl;
  cout << "Number of iterations: " << maxiter << endl;
  cout << "Tollerance achieved: " << tol << endl;
  if (verbose) {
    cout << "xk: " << endl;
    print_matrix(xk);
  }

  double relative_error = (x - xk).norm();
  cout << "Error norm (x - xk): " << relative_error << endl;
  if (verbose)
    print_matrix(x - xk);

  return 0;
}