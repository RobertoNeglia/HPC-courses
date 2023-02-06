#include "Heat.hpp"

// Main function.
int
main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int N      = 9;
  const unsigned int degree = 2;

  const double T      = 3.0;
  const double deltat = 0.0625;
  // const double theta  = 1.0; // implicit Euler (backward)
  // const double theta = 0.5; // Crank-Nicolson method (implicit)
  // E2.4 - no convergence, method is conditionally stable
  const double theta = 0.0; // explicit Euler (forward)

  Heat problem(N, degree, T, deltat, theta);

  problem.setup();
  problem.solve();

  return 0;
}