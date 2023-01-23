#ifndef POISSON_1D_HPP
#define POISSON_1D_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#define UNUSED(expr) \
    do {             \
      (void)(expr);  \
  } while (0)

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class Poisson2D {
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 2;

  // Diffusion coefficient.
  // In deal.ii, functions are implemented by deriving the dealii::Function
  // class, which provides an interface for the computation of function values
  // and their derivatives.
  class FunctionMu : public Function<dim> {
  public:
    // Constructor.
    FunctionMu() {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      UNUSED(p);
      return 1.0;
    }
  };

  // Reaction coefficient.
  class FunctionSigma : public Function<dim> {
  public:
    // Constructor
    FunctionSigma() {}

    // Evaluation
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      UNUSED(p);
      return 1.0;
    }
  };

  // Forcing term.
  class ForcingTerm : public Function<dim> {
  public:
    // Constructor.
    ForcingTerm() {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      return 1.0 - std::exp(p[0] + p[1]);
    }
  };

  // Neumann boundary condition constant.
  class FunctionQ : public Function<dim> {
  public:
    // Constructor.
    FunctionQ() {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      if (p[0] == 1 && (p[1] > 0 && p[1] < 1))
        return M_E * (std::exp(p[1]) - 1);
      if (p[1] == 1 && p[0] > 0 && p[0] < 1)
        return M_E * (std::exp(p[0]) - 1);
    }
  };

  // Exact solution
  class ExactSolution : public Function<dim> {
  public:
    // Constructor
    ExactSolution() {}

    // Evaluation
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      return (std::exp(p[0]) - 1) * (std::exp(p[1]) - 1);
    }

    // Evaluation of derivative
    virtual Tensor<1, dim>
    gradient(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      Tensor<1, dim> result;

      result[0] = std::exp(p[0]) * (std::exp(p[1]) - 1);
      result[1] = std::exp(p[1]) * (std::exp(p[0]) - 1);

      return result;
    }
  };

  // Constructor.
  Poisson2D(const unsigned int &N_, const unsigned int &r_) :
    N(N_), r(r_), mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
    mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), mesh(MPI_COMM_WORLD),
    pcout(std::cout, mpi_rank == 0) {}

  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve();

  // Output.
  void
  output() const;

  // Compute the error
  double
  compute_error(const VectorTools::NormType &norm_type) const;

protected:
  // N+1 is the number of elements.
  const unsigned int N;

  // Polynomial degree.
  const unsigned int r;

  // Number of MPI processes
  const unsigned int mpi_size;

  // This MPI process
  const unsigned int mpi_rank;

  // Diffusion coefficient.
  FunctionMu function_mu;

  // Robin boundary condition coefficient.
  FunctionSigma function_sigma;

  // Robin boundary condition constant.
  FunctionQ function_q;

  // Forcing term.
  ForcingTerm forcing_term;

  // Triangulation.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  // We use a unique_ptr here so that we can choose the type and degree of the
  // finite elements at runtime (the degree is a constructor parameter). The
  // class FiniteElement<dim> is an abstract class from which all types of
  // finite elements implemented by deal.ii inherit.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  // We use a unique_ptr here so that we can choose the type and order of the
  // quadrature formula at runtime (the order is a constructor parameter).
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // no sparsity pattern as a class member (its only needed in the setup method).

  // System matrix.
  TrilinosWrappers::SparseMatrix system_matrix;

  // System right-hand side.
  TrilinosWrappers::MPI::Vector system_rhs;

  // System solution.
  TrilinosWrappers::MPI::Vector solution;

  // Set of locally owned indices
  IndexSet locally_owned_dofs;

  ConditionalOStream pcout;

  ExactSolution exact_solution;
};

#endif