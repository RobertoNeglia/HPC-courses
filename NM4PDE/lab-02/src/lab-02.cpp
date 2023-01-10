#include <deal.II/base/convergence_table.h>

#include <iostream>

#include "Poisson1D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[]) {
  std::vector<unsigned int> N_values = {9, 19, 39, 79, 159};
  const unsigned int        r        = 2;

  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1" << std::endl;

  ConvergenceTable convergence_table;

    for (const unsigned int &N : N_values) {
      Poisson1D problem(N, r);

      problem.setup();
      problem.assemble();
      problem.solve();
      problem.output();

      // std::cout << "L2 error = " << problem.compute_error(VectorTools::L2_norm)
      //           << std::endl;
      // std::cout << "H1 error = " << problem.compute_error(VectorTools::H1_norm)
      //           << std::endl;

      const double h        = 1.0 / (N + 1);
      const double error_L2 = problem.compute_error(VectorTools::L2_norm);
      const double error_H1 = problem.compute_error(VectorTools::H1_norm);
      convergence_file << h << "," << error_L2 << "," << error_H1 << std::endl;

      convergence_table.add_value("h", h);
      convergence_table.add_value("eL2", error_L2);
      convergence_table.add_value("eH1", error_H1);
    }

  convergence_table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  convergence_table.set_scientific("eL2", true);
  convergence_table.set_scientific("eH1", true);

  convergence_table.write_text(std::cout);
  return 0;
}
