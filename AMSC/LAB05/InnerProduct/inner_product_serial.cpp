#define VERBOSE 0

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
inner_product(std::vector<double> a, std::vector<double> b) {
  if (a.size() != b.size())
    return 0;
  double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
      sum += a[i] * b[i];
    }

  return sum;
}

int
main(int argc, char **argv) {
  const int           N = 100003;
  double              result;
  std::vector<double> a(N), b(N);

    for (int i = 0; i < N; i++) {
      a[i] = i + 1;
      b[i] = N - i;
    }

    if (VERBOSE) {
      std::cout << "Vector a: " << std::endl;

        for (int i = 0; i < N; i++) {
          std::cout << a[i] << " - ";
        }

      std::cout << std::endl;

      std::cout << "Vector b: " << std::endl;

        for (int i = 0; i < N; i++) {
          std::cout << b[i] << " - ";
        }
      std::cout << std::endl;
  }

  const auto dt = timeit([&]() { result = inner_product(a, b); });

  std::cout << "Time elasped: " << dt << "[ms]" << std::endl;
  std::cout << "Result: " << result << std::endl;

  return 0;
}