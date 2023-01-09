#ifndef POISSON_1D_HPP
#define POISSON_1D_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
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

using namespace dealii;

class Poisson1D {
public:
  static constexpr unsigned int dim = 1;

  class DiffusionTerm : public Function<dim> {
  public:
    // Evaluate
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) {
      return 1.0;
    }

  protected:
  };

  class ForcingTerm : public Function<dim> {
  public:
    // Evaluate
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) {
      if (p[0] <= 1.0 / 8.0 || p[0] > 1.0 / 4.0)
        return 0.0;
      else
        return -1.0;
    }

  protected:
  };

  // Constructor.
  Poisson1D(const unsigned int &N_, const unsigned int &r_) : N(N_), r(r_) {}

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

  // Mesh
  Triangulation<dim> mesh;

  // Finite element space
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature rule
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler
  DoFHandler<dim> dof_handler;

  // Sparsity pattern of the matrix
  // part of the matrix that won't change if the matrix has to be updated
  SparsityPattern sparsity_pattern;

  // System matrix
  SparseMatrix<double> system_matrix;

  // System right-hand side
  Vector<double> system_rhs;

  // System solution
  Vector<double> system_solution;

  // Diffusion term
  DiffusionTerm diffusion_term;

  // Forcing term
  ForcingTerm forcing_term;
};

#endif