#include <omp.h>

#include <iostream>
#include <vector>

int
main() {
  const int           N = 12;
  std::vector<double> a(N), b(N), c(N);

    for (int i = 0; i < N; i++) {
      a[i] = 3;
      b[i] = 4;
      c[i] = 0;
    }

  int thread_id = 1, n_threads = 1, max_threads = 1;
#pragma omp parallel num_threads(N * 2) private(thread_id) shared(n_threads, max_threads)
  {
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif
#pragma omp for
      for (int i = 0; i < N; i++) {
#pragma omp critical
        std::cout << "Done by thread " << thread_id << std::endl;
        c[i] = a[i] + b[i];
      }

// #pragma omp critical
//     { std::cout << "Hi, from thread " << thread_id << std::endl; }
#ifdef _OPENMP
    n_threads   = omp_get_num_threads();
    max_threads = omp_get_max_threads();
#endif
  }

  std::cout << "Number of threads: " << n_threads << std::endl;
  std::cout << "Max number of threads: " << max_threads << std::endl;

  for (int i = 0; i < N; i++)
    std::cout << c[i] << std::endl;
}
