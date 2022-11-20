#include <mpi.h>

#include <chrono>
#include <functional>
#include <iostream>
#include <numeric>

auto
timeit(const std::function<void()> &f) {
  using namespace std::chrono;
  const auto start = high_resolution_clock::now();
  f();
  const auto end = high_resolution_clock::now();
  return duration_cast<milliseconds>(end - start).count();
}

double
parallel_inner_product(std::vector<double> a, std::vector<double> b) {
  if (a.size() != b.size())
    return 0;

  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // std::cout << "RANK: " << rank << ". HELLO!!\n Doing my little job..." << std::endl;

  // scatter the vectors into chunks
  auto N = a.size();
  MPI_Bcast(&N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

  std::vector<double> a_scatter, b_scatter;
  auto                scatter_size = N / size;

    // ASSUME N IS A MULTIPLE OF size (number of processors)
    if (!(N % size)) {
      if (!rank)
        std::cout << "N is divisible by size" << std::endl;
      a_scatter.resize(scatter_size);
      b_scatter.resize(scatter_size);

      MPI_Scatter(&a[0],
                  scatter_size,
                  MPI_DOUBLE,
                  &a_scatter[0],
                  scatter_size,
                  MPI_DOUBLE,
                  0,
                  MPI_COMM_WORLD);
      MPI_Scatter(&b[0],
                  scatter_size,
                  MPI_DOUBLE,
                  &b_scatter[0],
                  scatter_size,
                  MPI_DOUBLE,
                  0,
                  MPI_COMM_WORLD);

    } else { // N is not a multiple of size
      if (!rank)
        std::cout << "N is NOT divisible by size" << std::endl;
      scatter_size += (unsigned(rank) < N % size);
      a_scatter.resize(scatter_size);
      b_scatter.resize(scatter_size);

      std::vector<int> sendcounts(size), displs(size);

      displs[0] = 0;
        for (int i = 0; i < size; i++) {
          sendcounts[i] = N / size + (unsigned(i) < N % size);
            if (i > 0) {
              displs[i] = displs[i - 1] + sendcounts[i - 1];
          }
        }

      MPI_Scatterv(&a[0],
                   &sendcounts[0],
                   &displs[0],
                   MPI_DOUBLE,
                   &a_scatter[0],
                   scatter_size,
                   MPI_DOUBLE,
                   0,
                   MPI_COMM_WORLD);
      MPI_Scatterv(&b[0],
                   &sendcounts[0],
                   &displs[0],
                   MPI_DOUBLE,
                   &b_scatter[0],
                   scatter_size,
                   MPI_DOUBLE,
                   0,
                   MPI_COMM_WORLD);
    }

  // Now all processors have their chunk in a_scatter and b_scatter

  double partial_ip =
    std::inner_product(a_scatter.begin(), a_scatter.end(), b_scatter.begin(), 0.0);

    if (rank) {
      MPI_Send(&partial_ip, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
      double ip = partial_ip;
        for (int i = 1; i < size; i++) {
          MPI_Recv(&partial_ip, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ip += partial_ip;
        }

      return ip;
    }
  return 0;
}

int
main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int           N = 100003;
  std::vector<double> a(N), b(N);
  double              ip;

    for (int i = 0; i < N; i++) {
      a[i] = i + 1;
      b[i] = N - i;
    }

  // function call must be performed by all processors
  const auto dt = timeit([&]() { ip = parallel_inner_product(a, b); });

    if (!rank) {
      std::cout << "Time elapsed: " << dt << "[ms]" << std::endl;
      std::cout << "Inner product: " << ip << std::endl;
  }

  MPI_Finalize();
  return 0;
}
