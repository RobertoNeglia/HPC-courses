#ifndef POISSON_1D_HPP
#define POISSON_1D_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
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
  class DiffusionCoefficient : public Function<dim> {
  public:
    // Constructor.
    DiffusionCoefficient() {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      UNUSED(p);
      return 1.0;
    }
  };

  class ReactionCoefficient : public Function<dim> {
  public:
    // Constructor
    ReactionCoefficient() {}

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
      UNUSED(p);
      return 0.0;
    }
  };

  // Dirichlet boundary conditions.
  class FunctionG : public Function<dim> {
  public:
    // Constructor.
    FunctionG() {}

    // Evaluation.
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      UNUSED(p);
      return 1.0;
    }
  };

  // Neumann boundary conditions.
  class FunctionH : public Function<dim> {
  public:
    // Constructor.
    FunctionH() {}
    // Evaluation:
    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const {
      UNUSED(p);
      return 0.0;
    }
  };

  // Constructor.
  Poisson2D(const unsigned int &N_, const unsigned int &r_) : N(N_), r(r_) {}

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

protected:
  // N+1 is the number of elements.
  const unsigned int N;

  // Polynomial degree.
  const unsigned int r;

  // Diffusion coefficient.
  DiffusionCoefficient diffusion_coefficient;

  // Reaction coefficient
  ReactionCoefficient reaction_coefficient;

  // Forcing term.
  ForcingTerm forcing_term;

  // g(x).
  FunctionG function_g;

  // h(x).
  FunctionH function_h;

  // Triangulation.
  Triangulation<dim> mesh;

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

  // Sparsity pattern.
  SparsityPattern sparsity_pattern;

  // System matrix.
  SparseMatrix<double> system_matrix;

  // System right-hand side.
  Vector<double> system_rhs;

  // System solution.
  Vector<double> solution;
};

#endif