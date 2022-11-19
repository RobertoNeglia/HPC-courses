#include <mpi.h>

#include <chrono>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>

auto
timeit(const std::function<void()> &f) {
  using namespace std::chrono;
  const auto start = high_resolution_clock::now();
  f();
  const auto end = high_resolution_clock::now();
  return duration_cast<milliseconds>(end - start).count();
}

std::pair<double, double>
montecarlo(const std::function<double(double)> &f, unsigned long N) {
  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank) {
      MPI_Recv(&N, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        for (int i = 1; i < size; i++) {
          MPI_Send(&N, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
        }
    }

  // MPI_Bcast sends the message to all the processors, including itself. It's called by all
  // processors using the same arguments as the root, and on return the contents of root's buffer is
  // has been copied to all processors
  //   MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  // std::cout << "message broadcasted" << std::endl;

  // compute how many samples each processor has to evaluate
  const auto sample_block_size = N / size + (unsigned(rank) < (N % size));
  // std::cout << "block size: " << sample_block_size << std::endl;

  const auto                             seed = rank * rank * size * size;
  std::mt19937                           gen(seed);
  std::uniform_real_distribution<double> dist(-1., 1.);

  // std::cout << "engine instantiated, starting to generate samples..." << std::endl;

  double mean_f = 0., square_mean_f = 0.;
  double Q_n, var_Q_n;

    for (auto i = 0ul; i < sample_block_size; i++) {
      const auto fi = f(dist(gen));
      mean_f += fi;
      square_mean_f += fi * fi;
    }

  // std::cout << "summing up all partial sums..." << std::endl;

  // send partial sums to main proc
  //   double tmp1, tmp2;
  //   MPI_Reduce(&mean_f, &tmp1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //   mean_f = tmp1;
  //   MPI_Reduce(&square_mean_f, &tmp2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //   square_mean_f = tmp2;

    if (rank) {
      MPI_Send(&mean_f, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&square_mean_f, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    } else {
        for (int i = 1; i < size; i++) {
          double parial_mean, partial_square_mean;
          MPI_Recv(&parial_mean, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&partial_square_mean, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          mean_f += parial_mean;
          square_mean_f += partial_square_mean;
        }
    }

  // std::cout << "reduce done" << std::endl;

    if (!rank) {
      mean_f /= N;
      square_mean_f /= N;

      const double domain_measure = 2.0;
      Q_n                         = domain_measure * mean_f;
      var_Q_n = domain_measure * domain_measure / (N - 1) * (square_mean_f - mean_f * mean_f);
      return {Q_n, var_Q_n};
  } else
    return {0, 0};
}

int
main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const unsigned long N = 100000;

  const std::function<double(double)> f     = [](auto x) { return (x * x) - (4 * x) + 2; };
  const double                        exact = 14. / 3;
  const double                        tol   = 1.e-6;
  double                              var   = 1.;

    for (unsigned long N = 10; tol <= var; N *= 10) {
      std::pair<double, double> montecarlo_int;
      const auto                dt = timeit([&]() { montecarlo_int = montecarlo(f, N); });

        if (rank == 0) { // rank zero
          std::cout << "Exact result: " << exact << std::endl;
          std::cout << "Tolerance required: " << tol << std::endl;

          std::cout << "Number of samples: " << N << std::endl;
          std::cout << "Time elapsed: " << dt << "[ms]" << std::endl;
          std::cout << "Integration result: " << montecarlo_int.first << std::endl;
          std::cout << "Variance of integration result: " << montecarlo_int.second << std::endl;
          std::cout << "Error: " << std::abs(montecarlo_int.first - exact) << std::endl;
          var = montecarlo_int.second;
      }
    }
  MPI_Abort(MPI_COMM_WORLD, 0);

  MPI_Finalize();
  return 0;
}