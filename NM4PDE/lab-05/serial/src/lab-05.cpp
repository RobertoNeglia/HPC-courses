#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson3D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[]) {
  const unsigned int N = 39;
  const unsigned int r = 1;

  Poisson3D problem(N, r);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}