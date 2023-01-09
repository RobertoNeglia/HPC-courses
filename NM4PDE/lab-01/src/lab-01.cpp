#include <iostream>

#include "Poisson1D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[]) {
  unsigned int N = 39;
  unsigned int r = 1;
  Poisson1D    poisson_1D(N, r);

  poisson_1D.setup();
  poisson_1D.assemble();
  poisson_1D.solve();
  poisson_1D.output();

  return 0;
}
