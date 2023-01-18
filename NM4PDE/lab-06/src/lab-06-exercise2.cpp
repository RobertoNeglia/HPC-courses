#include "LinearElasticity.hpp"

// Main function.
int
main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int N      = 19;
  const unsigned int degree = 1;

  LinearElasticity problem(N, degree);

  problem.setup();
  problem.assemble_system();
  problem.solve_system();
  problem.output();

  return 0;
}