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

std::vector<double>
random_samples(unsigned long N) {
  std::random_device               rd; // will be used to obtain a seed for the random number engine
  std::mt19937                     gen(rd()); // standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-1.0, 1.0);
  std::vector<double>              samples(N);

    for (int i = 0; i < N; i++) {
      samples[i] = dis(gen);
    }
  return samples;
}

std::pair<double, double>
montecarlo(const std::function<double(double)> &f, unsigned long N) {
  const double domain_measure = 2.0; // volume of the domain
  double       mean_f, square_mean_f;
  double       Q_n, var_Q_n;

  // obtain samples in the range [-1.0, 1.0] to be evaluated
  const std::vector<double> samples = random_samples(N);

  //   mean_f =
  //     std::transform_reduce(samples.begin(), samples.end(), 0.0, std::plus<>(), [f](const double
  //     x) {
  //       return f(x);
  //     });

  //   square_mean_f =
  //     std::transform_reduce(samples.begin(), samples.end(), 0.0, std::plus<>(), [f](const double
  //     x) {
  //       double y = f(x);
  //       return y * y;
  //     });

  // perform the sum of all the evaluations of the samples
  mean_f        = 0.;
  square_mean_f = 0.;
    for (int i = 0; i < N; i++) {
      const double fi = f(samples[i]);
      mean_f += fi;
      square_mean_f += fi * fi;
    }

  mean_f /= N;
  square_mean_f /= N;

  // calculate the Monte Carlo approximation
  Q_n = domain_measure * mean_f;
  // calculate the variance
  var_Q_n = domain_measure * domain_measure / (N - 1) * (square_mean_f - mean_f * mean_f);

  return {Q_n, var_Q_n};
}

int
main(int argc, char **argv) {
  const std::function<double(double)> f     = [](auto x) { return (x * x) - (4 * x) + 2; };
  double                              var   = 1.0;
  const double                        tol   = 1.e-6;
  const double                        exact = 14. / 3;
  std::cout << "Exact result: " << exact << std::endl;

    for (int N = 10; var >= tol; N *= 10) {
      std::cout << "Number of samples: " << N << std::endl;
      std::pair<double, double> montecarlo_int;

      const auto dt = timeit([&]() { montecarlo_int = montecarlo(f, N); });

      std::cout << "Time elapsed: " << dt << "[ms]" << std::endl;
      std::cout << "Integration result: " << montecarlo_int.first << std::endl;
      std::cout << "Variance of integration result: " << montecarlo_int.second << std::endl;
      std::cout << "Error: " << std::abs(montecarlo_int.first - exact) << std::endl;
      var = montecarlo_int.second;
    }
  return 0;
}