#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "NonLinearDiffusion.hpp"

// Main function.
int
main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int N      = 39;
  const unsigned int degree = 1;

  NonLinearDiffusion problem(N, degree);

  problem.setup();
  problem.solve_newton();
  problem.output();

  return 0;
}