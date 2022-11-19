#include <mpi.h>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

// FIRST VERSION: ASSUME THE NUMBER OF THE PROCESSORS IS THE SAME AS THE ARRAY SIZE
template <size_t N>
std::array<double, N>
compute_gradient_v1(const std::function<double(const std::array<double, N> &)> &f,
                    const std::array<double, N>                                &y,
                    double                                                      h) {
  std::array<double, N> gradient;
  double                grad_comp;
  int                   rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<double, N> y_p(y);
  std::array<double, N> y_m(y);

  y_p[rank] += h;
  y_m[rank] -= h;

  grad_comp = (f(y_p) - f(y_m)) / (2 * h);

    if (rank) {
      MPI_Send(&grad_comp, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
      gradient[rank] = grad_comp;
        for (int i = 1; i < size; i++) {
          double rec_grad;
          MPI_Recv(&rec_grad, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          gradient[i] = rec_grad;
        }
    }

  return gradient;
}

// SECOND VERSION: THE ARRAY SIZE CAN BE ANY, COMPUTATION IS PERFORMED IN BLOCKS
template <size_t N>
std::array<double, N>
compute_gradient_v2(const std::function<double(const std::array<double, N> &)> &f,
                    const std::array<double, N>                                &y,
                    double                                                      h) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // calculate block size
  int sub_N = (N + size - 1) / size;

  std::array<double, N> gradient;
  std::vector<double>   sub_gradient(sub_N);

  // every processor computes his block
    for (int i = 0; i < sub_N; i++) {
      std::array<double, N> y_p(y);
      std::array<double, N> y_m(y);

      y_p[rank * sub_N + i] += h;
      y_m[rank * sub_N + i] -= h;

      sub_gradient[i] = (f(y_p) - f(y_m)) / (2 * h);
    }

    if (rank) {
      MPI_Send(&sub_gradient[0], sub_N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < sub_N; i++) {
          gradient[i] = sub_gradient[i];
        }
        for (int i = 1; i < size; i++) {
          std::vector<double> temp_gradient(sub_N);
          MPI_Recv(&temp_gradient[0], sub_N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < sub_N; j++) {
              if (i * sub_N + j < N)
                gradient[i * sub_N + j] = temp_gradient[j];
            }
        }
    }

  return gradient;
}

template <size_t N>
std::array<double, N>
compute_gradient_serial(const std::function<double(const std::array<double, N> &)> &f,
                        const std::array<double, N>                                &y,
                        double                                                      h) {
  std::array<double, N> gradient;

    for (int i = 0; i < N; i++) {
      std::array<double, N> y_p(y);
      std::array<double, N> y_m(y);
      y_p[i] += h;
      y_m[i] -= h;
      gradient[i] = (f(y_p) - f(y_m)) / (2 * h);
    }

  return gradient;
}

template <size_t N>
std::array<double, N>
diff_array(const std::array<double, N> &a, const std::array<double, N> &b) {
  std::array<double, N> diff;
  for (int i = 0; i < N; i++)
    diff[i] = std::abs(a[i] - b[i]);
  return diff;
}

template <size_t N>
bool
eq_rel_err(const std::array<double, N> &a, const std::array<double, N> &b, const double tol) {
    for (auto i : diff_array(a, b)) {
      if (i > tol)
        return false;
    }
  return true;
}

int
main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int                   rank, size;
  const int             n = 10000;
  std::array<double, n> x;

  for (int i = 0; i < n; i++)
    x[i] = (i + 1) * 2;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::function<double(const std::array<double, n> &)> f =
    [](auto x) { // just the sum of squares of every component
      return std::transform_reduce(x.begin(), x.end(), 0.0, std::plus<>(), [](const double x) {
        return x * x;
      });
    };

  std::array<double, n> gradient = compute_gradient_v2<n>(f, x, 1.e-6);

    if (!rank) {
      std::array<double, n> exact_gradient = compute_gradient_serial<n>(f, x, 1.e-6);

      std::cout << (eq_rel_err(exact_gradient, gradient, 1.e-6) ? "TEST PASSED" : "TEST NOT PASSED")
                << std::endl;
  }

  MPI_Finalize();
  return 0;
}