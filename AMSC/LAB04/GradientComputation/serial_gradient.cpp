#include <array>
#include <functional>
#include <iostream>
#include <numeric>

template <size_t N>
std::array<double, N>
compute_gradient(const std::function<double(const std::array<double, N> &)> &f,
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

int
main(int argc, char **argv) {
  const int             n = 2;
  std::array<double, n> x = {10, 1000000};

  const std::function<double(const std::array<double, n> &)> f =
    [](auto x) { // just the sum of squares of every component
      return x[0] + x[1];
    };

  double test = f(x);

  std::cout << "function evaluated in point [2,3]: " << test << std::endl;

  std::array<double, n> gradient = compute_gradient<n>(f, x, 1.e-6);

    for (auto d : gradient) {
      std::cout << d << std::endl;
    }

  return 0;
}